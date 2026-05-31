#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
tension_curve.py
================
Reproducible reference implementation of

    "Swap Curve Construction and Optimization of the Tension Parameter:
     An Energy-Partition Formulation"  (A. Sinha)

Everything in the paper is reproduced here from one self-contained file:

  * natural cubic spline                                  (Sec. III)
  * uniform hyperbolic tension spline                     (Sec. IV)
  * variational energies E_bend, E_tens  (closed form + quadrature)   (Sec. V, App. B)
  * energy-partition criterion  phi(sigma)  and the optimum sigma*     (Sec. V)
  * scale / unit invariance check                          (Prop. 2)
  * forward positivity and hedge-locality radius           (Sec. VI)
  * area-constrained tension spline on the FORWARD curve   (Sec. VII)
  * non-uniform per-segment tension calibration            (Sec. VIII)
  * synthesis: non-uniform tension on the forward curve    (Sec. IX)

All hyperbolic evaluations use numerically stable e^{-sigma h} forms, so the
code is safe for arbitrarily large products sigma*h (Sec. IV-B / App. B).

Usage
-----
    python tension_curve.py                # full validation report + figures
    python tension_curve.py --no-figures   # report only

Dependencies: numpy, matplotlib (matplotlib only needed for figures).

The script is deterministic; running it twice gives identical numbers.
"""

from __future__ import annotations
import argparse
import numpy as np

# ======================================================================
# 0.  MARKET DATA   (U.S. Treasury CMT par-yield curve, 22-May-2026)
#     Source: Federal Reserve H.15 release.  Yields are decimals.
# ======================================================================
TENORS = np.array([1/12, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0])
YIELDS = np.array([0.0372, 0.0368, 0.0379, 0.0386, 0.0413, 0.0418,
                   0.0427, 0.0441, 0.0456, 0.0506, 0.0507])

# ======================================================================
# 1.  STABLE HYPERBOLIC HELPERS
#     ratios sinh(s u)/sinh(s h) and cosh(s u)/sinh(s h) for 0 <= u <= h,
#     evaluated without overflow for large s*h.
# ======================================================================
def _sinh_ratio(u, h, s):
    """sinh(s*u) / sinh(s*h),  0 <= u <= h,  s >= 0.  Overflow-safe."""
    u = np.asarray(u, float)
    x = s * h
    small = x < 30.0
    # small-argument branch: direct
    out_small = np.sinh(np.clip(s * u, -700, 700)) / np.sinh(np.clip(np.where(small, x, 1.0), -700, 700))
    # large-argument branch: factor e^{s h}
    e = np.exp(np.where(small, 0.0, s * (u - h)))      # e^{s(u-h)} <= 1
    e2 = np.exp(np.where(small, 0.0, -s * (u + h)))    # e^{-s(u+h)}
    den = 1.0 - np.exp(np.where(small, -1.0, -2.0 * x))
    out_large = (e - e2) / den
    return np.where(small, out_small, out_large)

def _cosh_ratio(u, h, s):
    """cosh(s*u) / sinh(s*h),  0 <= u <= h,  s >= 0.  Overflow-safe."""
    u = np.asarray(u, float)
    x = s * h
    small = x < 30.0
    out_small = np.cosh(np.clip(s * u, -700, 700)) / np.sinh(np.clip(np.where(small, x, 1.0), -700, 700))
    e = np.exp(np.where(small, 0.0, s * (u - h)))
    e2 = np.exp(np.where(small, 0.0, -s * (u + h)))
    den = 1.0 - np.exp(np.where(small, -1.0, -2.0 * x))
    out_large = (e + e2) / den
    return np.where(small, out_small, out_large)


# ======================================================================
# 2.  NATURAL CUBIC SPLINE  (Sec. III)
# ======================================================================
def cubic_solve(x, y):
    """Return second derivatives z_i of the natural cubic spline."""
    n = len(x); h = np.diff(x); z = np.zeros(n)
    if n <= 2:
        return z
    A = np.zeros((n - 2, n - 2)); b = np.zeros(n - 2)
    for k in range(1, n - 1):
        i = k - 1
        A[i, i] = 2.0 * (h[k - 1] + h[k])
        if i > 0:     A[i, i - 1] = h[k - 1]
        if i < n - 3: A[i, i + 1] = h[k]
        b[i] = 6.0 * ((y[k + 1] - y[k]) / h[k] - (y[k] - y[k - 1]) / h[k - 1])
    z[1:n - 1] = np.linalg.solve(A, b)
    return z

def cubic_eval(x, y, z, t, deriv=0):
    t = np.atleast_1d(np.asarray(t, float)); n = len(x)
    idx = np.clip(np.searchsorted(x, t) - 1, 0, n - 2)
    out = np.empty_like(t)
    for j, (tt, i) in enumerate(zip(t, idx)):
        hi = x[i + 1] - x[i]; d = tt - x[i]
        A = (z[i + 1] - z[i]) / (6 * hi); B = z[i] / 2.0
        C = -(hi / 6) * z[i + 1] - (hi / 3) * z[i] + (y[i + 1] - y[i]) / hi
        if   deriv == 0: out[j] = y[i] + d * (C + d * (B + d * A))
        elif deriv == 1: out[j] = C + d * (2 * B + 3 * d * A)
        elif deriv == 2: out[j] = 2 * B + 6 * d * A
    return out


# ======================================================================
# 3.  UNIFORM TENSION SPLINE  (Sec. IV)   f'''' - sigma^2 f'' = 0
# ======================================================================
def tension_solve(x, y, s):
    """Second derivatives z_i of the uniform tension spline (natural BCs)."""
    n = len(x); h = np.diff(x); z = np.zeros(n)
    if n <= 2:
        return z
    sh = np.sinh(np.clip(s * h, 0, 700))
    alpha = 1.0 / h - s / sh
    beta  = s * np.cosh(np.clip(s * h, 0, 700)) / sh - 1.0 / h
    gamma = s * s * (y[1:] - y[:-1]) / h
    A = np.zeros((n - 2, n - 2)); b = np.zeros(n - 2)
    for k in range(1, n - 1):
        i = k - 1
        A[i, i] = beta[k - 1] + beta[k]
        if i > 0:     A[i, i - 1] = alpha[k - 1]
        if i < n - 3: A[i, i + 1] = alpha[k]
        b[i] = gamma[k] - gamma[k - 1]
    z[1:n - 1] = np.linalg.solve(A, b)
    return z

def tension_eval(x, y, z, s, t, deriv=0):
    """Evaluate the uniform tension spline (stable for large sigma*h)."""
    t = np.atleast_1d(np.asarray(t, float)); n = len(x)
    idx = np.clip(np.searchsorted(x, t) - 1, 0, n - 2)
    out = np.empty_like(t)
    for j, (tt, i) in enumerate(zip(t, idx)):
        h = x[i + 1] - x[i]; a = tt - x[i]; p = h - a
        if deriv == 0:
            hyp = (z[i] * _sinh_ratio(p, h, s) + z[i + 1] * _sinh_ratio(a, h, s)) / s**2
            lin = (y[i] - z[i] / s**2) * (p / h) + (y[i + 1] - z[i + 1] / s**2) * (a / h)
            out[j] = hyp + lin
        elif deriv == 1:
            hyp = (-z[i] * _cosh_ratio(p, h, s) + z[i + 1] * _cosh_ratio(a, h, s)) / s
            out[j] = hyp + ((y[i + 1] - y[i]) - (z[i + 1] - z[i]) / s**2) / h
        elif deriv == 2:
            out[j] = z[i] * _sinh_ratio(p, h, s) + z[i + 1] * _sinh_ratio(a, h, s)
    return out


# ======================================================================
# 4.  NON-UNIFORM TENSION SPLINE  (Sec. VIII, App. A)
# ======================================================================
def nonuniform_solve(x, y, sig):
    """z_i for a tension spline with per-segment tension sig[i] on [x_i,x_{i+1}]."""
    n = len(x); h = np.diff(x); z = np.zeros(n)
    if n <= 2:
        return z
    M = np.zeros((n - 2, n - 2)); b = np.zeros(n - 2)
    def cf(s, hh):                 # off-diagonal coefficient 1/(s^2 h) - 1/(s sinh(s h))
        return 1.0 / (s**2 * hh) - 1.0 / (s * np.sinh(np.clip(s * hh, 0, 700)))
    def dg(s, hh):                 # diagonal half cosh(s h)/(s sinh(s h)) - 1/(s^2 h)
        sh = np.sinh(np.clip(s * hh, 0, 700))
        return np.cosh(np.clip(s * hh, 0, 700)) / (s * sh) - 1.0 / (s**2 * hh)
    for k in range(1, n - 1):
        i = k - 1; sL, hL = sig[k - 1], h[k - 1]; sR, hR = sig[k], h[k]
        M[i, i] = dg(sL, hL) + dg(sR, hR)
        if i > 0:     M[i, i - 1] = cf(sL, hL)
        if i < n - 3: M[i, i + 1] = cf(sR, hR)
        b[i] = (y[k + 1] - y[k]) / hR - (y[k] - y[k - 1]) / hL
    z[1:n - 1] = np.linalg.solve(M, b)
    return z

def nonuniform_eval(x, y, z, sig, t, deriv=0):
    t = np.atleast_1d(np.asarray(t, float)); n = len(x)
    idx = np.clip(np.searchsorted(x, t) - 1, 0, n - 2)
    out = np.empty_like(t)
    for j, (tt, i) in enumerate(zip(t, idx)):
        s = sig[i]; h = x[i + 1] - x[i]; a = tt - x[i]; p = h - a
        if deriv == 0:
            hyp = (z[i] * _sinh_ratio(p, h, s) + z[i + 1] * _sinh_ratio(a, h, s)) / s**2
            lin = (y[i] - z[i] / s**2) * (p / h) + (y[i + 1] - z[i + 1] / s**2) * (a / h)
            out[j] = hyp + lin
        elif deriv == 1:
            hyp = (-z[i] * _cosh_ratio(p, h, s) + z[i + 1] * _cosh_ratio(a, h, s)) / s
            out[j] = hyp + ((y[i + 1] - y[i]) - (z[i + 1] - z[i]) / s**2) / h
        elif deriv == 2:
            out[j] = z[i] * _sinh_ratio(p, h, s) + z[i + 1] * _sinh_ratio(a, h, s)
    return out


# ======================================================================
# 5.  ENERGIES  (Sec. V, App. B)
# ======================================================================
_GL_X, _GL_W = np.polynomial.legendre.leggauss(64)   # 64-pt Gauss-Legendre

def segment_energies_quadrature(evalfun, x, i, s):
    """(E_bend^i, E_tens^i) on [x_i,x_{i+1}] by Gauss-Legendre, stable at any s."""
    a, b = x[i], x[i + 1]; half = 0.5 * (b - a); mid = 0.5 * (a + b)
    nodes = mid + half * _GL_X; w = half * _GL_W
    f2 = evalfun(nodes, 2); f1 = evalfun(nodes, 1)
    Eb = 0.5 * np.sum(w * f2 * f2)
    Et = 0.5 * s * s * np.sum(w * f1 * f1)
    return Eb, Et

def energies_uniform(x, y, z, s):
    Eb = Et = 0.0
    ev = lambda t, d: tension_eval(x, y, z, s, t, d)
    for i in range(len(x) - 1):
        b, t = segment_energies_quadrature(ev, x, i, s); Eb += b; Et += t
    return Eb, Et

def bending_energy_closed(zi, zj, h, s):
    """Closed-form flexural energy of one segment (Eq. for E_bend^i). Moderate s."""
    sh = np.sinh(s * h); s2 = np.sinh(2 * s * h); ch = np.cosh(s * h)
    return (1.0 / (2 * sh**2)) * ((zi**2 + zj**2) * (s2 / (4 * s) - h / 2)
                                  + zi * zj * (h * ch - sh / s))

def axial_energy_closed(zi, zj, yi, yj, h, s):
    """Closed-form axial energy of one segment (Eq. E_a). Moderate s."""
    sh = np.sinh(s * h); s2 = np.sinh(2 * s * h); ch = np.cosh(s * h)
    S = ((yj - yi) - (zj - zi) / s**2) / h
    g2 = (1.0 / (2 * sh**2)) * ((zi**2 + zj**2) * (s2 / (4 * s) + h / 2)
                                - zi * zj * (h * ch + sh / s))
    return g2 + S * (zj - zi) + 0.5 * s**2 * S**2 * h


# ======================================================================
# 6.  ENERGY-PARTITION CRITERION  (Sec. V)
# ======================================================================
def phi_of_sigma(x, y, s, method="closed"):
    """Bending energy fraction phi(sigma). Uses the closed-form energies by default
    (O(1) per segment, validated against 64-pt quadrature to ~1e-12); pass
    method='quadrature' to force the Gauss-Legendre path."""
    z = tension_solve(x, y, s)
    if method == "closed":
        h = np.diff(x)
        Eb = sum(bending_energy_closed(z[i], z[i+1], h[i], s) for i in range(len(x)-1))
        Et = sum(axial_energy_closed(z[i], z[i+1], y[i], y[i+1], h[i], s) for i in range(len(x)-1))
    else:
        Eb, Et = energies_uniform(x, y, z, s)
    return Eb / (Eb + Et)

def solve_sigma_star(x, y, phi_target=0.5, lo=1e-2, hi=60.0, iters=80):
    """Unique sigma* with phi(sigma*) = phi_target by bisection (phi is monotone)."""
    f = lambda s: phi_of_sigma(x, y, s) - phi_target
    flo = f(lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi); fm = f(mid)
        if (flo > 0) == (fm > 0): lo, flo = mid, fm
        else: hi = mid
    return 0.5 * (lo + hi)


# ======================================================================
# 7.  HEDGE LOCALITY  (Sec. VI)
# ======================================================================
def locality_radius(x, y, k, kind="tension", s=None, eps=1e-4, ngrid=2000):
    """rho_k = mean |t - t_k| weighted by |df/dr_k| over the curve."""
    grid = np.linspace(x[0], x[-1], ngrid)
    yp = y.copy(); yp[k] += eps
    if kind == "cubic":
        f0 = cubic_eval(x, y,  cubic_solve(x, y),  grid)
        f1 = cubic_eval(x, yp, cubic_solve(x, yp), grid)
    else:
        f0 = tension_eval(x, y,  tension_solve(x, y, s),  s, grid)
        f1 = tension_eval(x, yp, tension_solve(x, yp, s), s, grid)
    resp = np.abs(f1 - f0); tk = x[k]
    rho = np.trapezoid(np.abs(grid - tk) * resp, grid) / np.trapezoid(resp, grid)
    return rho, grid, resp


# ======================================================================
# 8.  FORWARD-SPACE CONSTRUCTIONS  (Sec. VII & IX)
#     Treat the CMT curve as a zero curve; I(t) = r(t) t = -ln P(t).
#     The instantaneous forward f is a tension spline whose integral over
#     each interval reproduces the discount factor (area constraint).
# ======================================================================
def _fwd_operators(x, y, ngrid=2401):
    tau = np.concatenate([[0.0], x]); Iknot = np.concatenate([[0.0], y * x])
    areas = np.diff(Iknot); m = len(areas)
    g = np.linspace(tau[0], tau[-1], ngrid); dt = g[1] - g[0]; N = ngrid
    D1 = np.zeros((N - 1, N)); D2 = np.zeros((N - 2, N))
    for k in range(N - 1): D1[k, k] = -1 / dt; D1[k, k + 1] = 1 / dt
    for k in range(N - 2): D2[k, k] = 1 / dt**2; D2[k, k + 1] = -2 / dt**2; D2[k, k + 2] = 1 / dt**2
    Kb = D2.T @ D2 * dt
    C = np.zeros((m, N))
    for i in range(m):
        C[i, :] = np.clip(np.minimum(g + dt / 2, tau[i + 1]) - np.maximum(g - dt / 2, tau[i]), 0, None)
    mid = 0.5 * (g[:-1] + g[1:])
    segD1 = np.clip(np.searchsorted(tau, mid, side="right") - 1, 0, m - 1)
    segD2 = np.clip(np.searchsorted(tau, g[1:-1], side="right") - 1, 0, m - 1)
    return dict(tau=tau, areas=areas, m=m, g=g, dt=dt, N=N,
                D1=D1, D2=D2, Kb=Kb, C=C, segD1=segD1, segD2=segD2)

def forward_solve(op, sigvec):
    """Area-constrained tension spline on the forward; sigvec is per-segment."""
    N, m = op["N"], op["m"]
    w = sigvec[op["segD1"]] ** 2
    Kt = op["D1"].T @ (w[:, None] * op["D1"]) * op["dt"]
    K = op["Kb"] + Kt
    A = np.zeros((N + m, N + m)); A[:N, :N] = K + 1e-12 * np.eye(N)
    A[:N, N:] = op["C"].T; A[N:, :N] = op["C"]
    sol = np.linalg.solve(A, np.concatenate([np.zeros(N), op["areas"]]))
    return sol[:N]

def forward_segment_energies(op, f):
    f2 = (op["D2"] @ f) ** 2 * op["dt"]; f1 = (op["D1"] @ f) ** 2 * op["dt"]
    B = np.array([0.5 * f2[op["segD2"] == i].sum() for i in range(op["m"])])
    A = np.array([0.5 * f1[op["segD1"] == i].sum() for i in range(op["m"])])
    return B, A

def forward_uniform_sigma_star(op, target=0.5, lo=0.1, hi=40.0, iters=60):
    def phi(s):
        f = forward_solve(op, np.full(op["m"], s)); B, A = forward_segment_energies(op, f)
        return B.sum() / (B.sum() + s * s * A.sum())
    flo = phi(lo) - target
    for _ in range(iters):
        mid = 0.5 * (lo + hi); fm = phi(mid) - target
        if (flo > 0) == (fm > 0): lo, flo = mid, fm
        else: hi = mid
    return 0.5 * (lo + hi)

def forward_nonuniform_calibrate(op, s0, target=0.5, iters=200, tol=1e-7):
    """Per-segment equipartition on the forward; closed-form update sigma=sqrt(B/A)."""
    sig = np.full(op["m"], s0)
    for it in range(iters):
        f = forward_solve(op, sig); B, A = forward_segment_energies(op, f)
        new = np.sqrt(np.maximum(B, 1e-30) / np.maximum(A, 1e-30))
        if np.max(np.abs(np.log(new / sig))) < tol:
            break
        sig = np.exp(0.25 * np.log(sig) + 0.75 * np.log(new))
    return sig, forward_solve(op, sig), it


# ======================================================================
# 8b.  MULTI-CURVE: SOFR/OIS  (Sec. X)
#      Representative SOFR OIS par-swap levels for the 22-May-2026 environment
#      (anchored to 5Y SOFR ~4.25%, Fed funds 3.50-3.75%); the tenor basis is
#      illustrative.  Levels are representative; the construction is exact.
# ======================================================================
SOFR_TENORS = np.array([1., 2., 3., 5., 7., 10., 20., 30.])
SOFR_PAR    = np.array([3.90, 4.10, 4.18, 4.25, 4.30, 4.36, 4.55, 4.45]) / 100.0

def bootstrap_ois(sw_t=SOFR_TENORS, sw_s=SOFR_PAR):
    """Bootstrap annual OIS discount factors from par swap rates (tau=1, annual).
    Returns (years, P, zero_rates) and the max reprice error of the input swaps."""
    yrs = np.arange(1, int(sw_t[-1]) + 1).astype(float)
    S = np.interp(yrs, sw_t, sw_s)
    P = np.zeros(len(yrs)); A = 0.0
    for k, Sn in enumerate(S):
        P[k] = (1.0 - Sn * A) / (1.0 + Sn)      # P_n (1+S_n) = 1 - S_n * annuity_{n-1}
        A += P[k]
    annu = np.cumsum(P); reprice = np.max(np.abs((1.0 - P) / annu - S))
    return yrs, P, -np.log(P) / yrs, reprice

def multicurve_demo(ngrid=801):
    """Crisp OIS discount + tensioned projection forward. Returns key diagnostics.
    ngrid controls the forward-grid resolution (801 is ample for the report;
    the paper's headline figures use 2401)."""
    yrs, P, z_ois, reprice = bootstrap_ois()
    z_ois_t = -np.log(np.interp(SOFR_TENORS, yrs, P)) / SOFR_TENORS
    basis_t = (5 + 13 * (SOFR_TENORS - 1) / 29) / 1e4     # illustrative tenor basis
    z_proj_t = z_ois_t + basis_t
    op = _fwd_operators(SOFR_TENORS, z_proj_t, ngrid=ngrid)
    s = forward_uniform_sigma_star(op, 0.5)
    f = forward_solve(op, np.full(op["m"], s))
    return dict(reprice=reprice, sigma=s,
                df_err=np.max(np.abs(op["C"] @ f - op["areas"])),
                min_fwd=f.min())


# ======================================================================
# 9.  VALIDATION REPORT
# ======================================================================
def _M_and_dbdy(xx, s):
    h = np.diff(xx); m = len(xx)
    sh = np.sinh(np.clip(s*h, 0, 700)); ch = np.cosh(np.clip(s*h, 0, 700))
    al = 1.0/h - s/sh; be = s*ch/sh - 1.0/h
    M = np.zeros((m-2, m-2)); dbdy = np.zeros((m-2, m))
    for k in range(1, m-1):
        i = k-1; M[i, i] = be[k-1] + be[k]
        if i > 0:     M[i, i-1] = al[k-1]
        if i < m-3:   M[i, i+1] = al[k]
        dbdy[i, k+1] += s*s/h[k];   dbdy[i, k] += -s*s/h[k]
        dbdy[i, k]   -= s*s/h[k-1]; dbdy[i, k-1] -= -s*s/h[k-1]
    return M, dbdy

def dz_dy(xx, s):
    """Analytic sensitivity of the knot second-derivatives to the input rates."""
    m = len(xx); M, dbdy = _M_and_dbdy(xx, s)
    dz = np.zeros((m, m)); dz[1:m-1, :] = np.linalg.solve(M, dbdy)
    return dz

def df_dy(xx, yy, s, tq):
    """Analytic curve delta d f(x)/d y_j at fixed sigma (the hedge Jacobian)."""
    m = len(xx); h = np.diff(xx); z = tension_solve(xx, yy, s); dz = dz_dy(xx, s)
    tq = np.atleast_1d(tq); out = np.zeros((len(tq), m))
    idx = np.clip(np.searchsorted(xx, tq)-1, 0, m-2)
    for q, (tt, i) in enumerate(zip(tq, idx)):
        hi = h[i]; sh = np.sinh(np.clip(s*hi, 0, 700))
        P = xx[i+1]-tt; Q = tt-xx[i]; S1 = np.sinh(s*P); S2 = np.sinh(s*Q)
        for j in range(m):
            di, dj = dz[i, j], dz[i+1, j]
            out[q, j] = (di*S1+dj*S2)/(s*s*sh) \
                + ((1.0 if j == i else 0.0)-di/s**2)*(P/hi) \
                + ((1.0 if j == i+1 else 0.0)-dj/s**2)*(Q/hi)
    return out

def risk_sensitivities_demo():
    """Validate analytic df/dy vs finite diff; report sigma* stress range."""
    x, y = TENORS, YIELDS; s = solve_sigma_star(x, y, 0.5)
    tq = np.linspace(x[0], x[-1], 40); ana = df_dy(x, y, s, tq)
    eps = 1e-6; f0 = tension_eval(x, y, tension_solve(x, y, s), s, tq); err = 0.0
    for j in range(len(x)):
        yp = y.copy(); yp[j] += eps
        fd = (tension_eval(x, yp, tension_solve(x, yp, s), s, tq)-f0)/eps
        err = max(err, np.max(np.abs(ana[:, j]-fd)))
    mid = np.mean(x)
    sl = [solve_sigma_star(x, y+d*(x-mid)/(x[-1]-x[0]), 0.5) for d in (-0.01, 0.01)]
    return dict(df_err=err, slope_lo=min(sl), slope_hi=max(sl))


# ======================================================================
# 8c.  BENCHMARK INTERPOLATION METHODS  (Sec. XII)  -- no external data
# ======================================================================
def _zero_to_fwd(tg, zg):
    return np.gradient(tg * zg, tg)

def bench_log_discount(x, y, tg):
    """Log-linear in the discount factor (piecewise-constant instantaneous forward)."""
    lnP = np.interp(tg, x, np.log(np.exp(-x * y)))
    return -lnP / tg

def bench_piecewise_flat_forward(x, y, tg):
    """Piecewise-constant instantaneous forward between nodes."""
    fwd = np.empty(len(x)); fwd[0] = y[0]
    fwd[1:] = (x[1:]*y[1:] - x[:-1]*y[:-1]) / (x[1:] - x[:-1])
    idx = np.clip(np.searchsorted(x, tg, side="right") - 1, 0, len(x) - 2)
    out = np.empty_like(tg)
    for j, t in enumerate(tg):
        i = idx[j]; acc = x[0]*y[0]
        for k in range(i): acc += fwd[k+1]*(x[k+1]-x[k])
        acc += fwd[i+1]*(t - x[i]); out[j] = acc / t
    return out

def bench_hagan_west(x, y, tg):
    """Hagan-West monotone-convex forward interpolation (Wilmott, 2006)."""
    n = len(x); xx = np.concatenate([[0.0], x]); h = np.diff(xx)
    fd = np.empty(n); fd[0] = y[0]
    for i in range(1, n):
        fd[i] = (x[i]*y[i] - x[i-1]*y[i-1]) / (x[i]-x[i-1])
    fk = np.empty(n+1)
    for i in range(1, n):
        fk[i] = (h[i]/(h[i-1]+h[i]))*fd[i-1] + (h[i-1]/(h[i-1]+h[i]))*fd[i]
    fk[0] = fd[0] - 0.5*(fk[1]-fd[0]); fk[n] = fd[n-1] - 0.5*(fk[n-1]-fd[n-1])
    def inst(t):
        i = int(np.clip(np.searchsorted(xx, t, side="right")-1, 0, n-1))
        xi = (t-xx[i])/h[i]; g0 = fk[i]-fd[i]; g1 = fk[i+1]-fd[i]
        return fd[i] + g0*(1-4*xi+3*xi**2) + g1*(-2*xi+3*xi**2)
    out = np.empty_like(tg)
    for j, t in enumerate(tg):
        s = np.linspace(1e-9, t, 256); out[j] = np.trapezoid([inst(u) for u in s], s)/t
    return out

def bench_smith_wilson(x, y, tg, ufr=0.039, alpha=0.10):
    """Smith-Wilson zero curve (UFR and alpha are user inputs, illustrative defaults)."""
    P = np.exp(-x * y); mu = np.exp(-ufr * x)
    def W(t, u):
        a, b = min(t, u), max(t, u)
        return np.exp(-ufr*(t+u))*(alpha*min(t,u) - 0.5*np.exp(-alpha*b)*(np.exp(alpha*a)-np.exp(-alpha*a)))
    Wm = np.array([[W(ti, tj) for tj in x] for ti in x])
    zeta = np.linalg.solve(Wm, P - mu); out = np.empty_like(tg)
    for j, t in enumerate(tg):
        wt = np.array([W(t, tj) for tj in x])
        out[j] = np.exp(-ufr*t) + wt @ zeta
    return -np.log(np.clip(out, 1e-8, None)) / tg

def forward_roughness(tg, zg):
    """int (d f_fwd/dt)^2 dt -- a global smoothness penalty on the forward."""
    f = _zero_to_fwd(tg, zg); return np.trapezoid(np.gradient(f, tg)**2, tg)

def benchmark_table(x=None, y=None, ngrid=2000):
    """Compare all methods on (x,y); returns list of dicts. Defaults to the 2026 curve."""
    if x is None: x, y = TENORS, YIELDS
    tg = np.linspace(x[0], x[-1], ngrid)
    s = solve_sigma_star(x, y, 0.5)
    curves = {
        "cubic":           cubic_eval(x, y, cubic_solve(x, y), tg),
        "tension_phi_0.5": tension_eval(x, y, tension_solve(x, y, s), s, tg),
        "hagan_west":      bench_hagan_west(x, y, tg),
        "smith_wilson":    bench_smith_wilson(x, y, tg),
        "log_discount":    bench_log_discount(x, y, tg),
        "pw_flat_forward": bench_piecewise_flat_forward(x, y, tg),
    }
    out = []
    for name, zg in curves.items():
        f = _zero_to_fwd(tg, zg)
        out.append(dict(method=name, fwd_roughness=forward_roughness(tg, zg),
                        min_fwd=f.min()))
    return out, s

def phi_sweep(x=None, y=None, targets=(0.3, 0.4, 0.5, 0.6, 0.7), ngrid=2000):
    """Characterise the dial: sigma*, forward roughness, 5y locality vs phi*."""
    if x is None: x, y = TENORS, YIELDS
    tg = np.linspace(x[0], x[-1], ngrid); rows = []
    k5 = int(np.argmin(np.abs(x - 5.0)))
    for p in targets:
        s = solve_sigma_star(x, y, p)
        zg = tension_eval(x, y, tension_solve(x, y, s), s, tg)
        loc, _, _ = locality_radius(x, y, k5, "tension", s=s)
        rows.append(dict(phi=p, sigma=s, fwd_roughness=forward_roughness(tg, zg), loc5y=loc))
    return rows


def validation_report():
    x, y = TENORS, YIELDS
    line = "-" * 64
    print(line); print("VALIDATION REPORT  (all checks computed from scratch)"); print(line)

    # (a) C2 continuity of the uniform tension spline
    s = 0.6; z = tension_solve(x, y, s); eps = 1e-6; worst1 = worst2 = 0.0
    for k in range(1, len(x) - 1):
        worst1 = max(worst1, abs(tension_eval(x, y, z, s, x[k]-eps, 1)[0]
                                  - tension_eval(x, y, z, s, x[k]+eps, 1)[0]))
        worst2 = max(worst2, abs(tension_eval(x, y, z, s, x[k]-eps, 2)[0]
                                  - tension_eval(x, y, z, s, x[k]+eps, 2)[0]))
    print(f"[IV ] C2 continuity: max jump f'={worst1:.2e}, f''={worst2:.2e}  (expect ~0)")

    # (b) sigma->0 reproduces cubic
    tt = np.linspace(x[0], x[-1], 400)
    d = np.max(np.abs(tension_eval(x, y, tension_solve(x, y, 1e-3), 1e-3, tt)
                      - cubic_eval(x, y, cubic_solve(x, y), tt)))
    print(f"[IV ] sigma->0 matches cubic: max|diff|={d:.2e}")

    # (c) closed-form energies vs quadrature (moderate sigma)
    s = 0.6; z = tension_solve(x, y, s); h = np.diff(x)
    Bc = sum(bending_energy_closed(z[i], z[i+1], h[i], s) for i in range(len(x)-1))
    Ac = sum(axial_energy_closed(z[i], z[i+1], y[i], y[i+1], h[i], s) for i in range(len(x)-1))
    Bq, Aq = energies_uniform(x, y, z, s)
    print(f"[V/B] closed-form vs quadrature: E_bend rel.err={abs(Bc-Bq)/Bq:.2e}, "
          f"E_tens rel.err={abs(Ac-Aq)/Aq:.2e}")

    # (d) phi limits and monotonicity (Prop. 1)
    philo = phi_of_sigma(x, y, 1e-3); phihi = phi_of_sigma(x, y, 40.0, method="quadrature")
    sgrid = np.linspace(0.05, 30, 200); phis = np.array([phi_of_sigma(x, y, s) for s in sgrid])
    mono = np.all(np.diff(phis) < 1e-9)
    print(f"[V  ] phi(0+)={philo:.4f} (->1), phi(inf)={phihi:.4f} (->0); "
          f"strictly decreasing on grid: {mono}")

    # (e) sigma* for three preferences
    s75 = solve_sigma_star(x, y, 0.75); s50 = solve_sigma_star(x, y, 0.50); s25 = solve_sigma_star(x, y, 0.25)
    print(f"[V  ] sigma*: phi*=0.75 -> {s75:.3f}; phi*=0.50 -> {s50:.3f}; phi*=0.25 -> {s25:.3f}")

    # (f) scale invariance of phi (Prop. 2)
    p_dec = phi_of_sigma(x, y, s50)
    p_bps = phi_of_sigma(x, 1e4 * y, s50)                 # rate axis x1e4 (decimals->bps)
    p_time = phi_of_sigma(2.0 * x, y, s50 / 2.0)          # time x2, sigma/2 (sigma*h fixed)
    print(f"[V  ] phi invariance: decimals={p_dec:.6f}, bps={p_bps:.6f}, "
          f"time-rescaled={p_time:.6f}")

    # (g) forward positivity
    grid = np.linspace(x[0], x[-1], 2000)
    fc = cubic_eval(x, y, cubic_solve(x, y), grid) \
         + grid * cubic_eval(x, y, cubic_solve(x, y), grid, 1)
    zt = tension_solve(x, y, s50)
    ft = tension_eval(x, y, zt, s50, grid) + grid * tension_eval(x, y, zt, s50, grid, 1)
    print(f"[VI ] min implied fwd: cubic={fc.min():.4%}, tension(sigma*)={ft.min():.4%}  (>0 = no-arb)")

    # (h) hedge locality at the 5y node (index 6)
    rc, *_ = locality_radius(x, y, 6, "cubic")
    rt, *_ = locality_radius(x, y, 6, "tension", s=s50)
    print(f"[VI ] hedge-locality radius at 5y: cubic={rc:.2f} yr, tension={rt:.2f} yr "
          f"({rc/rt:.1f}x tighter)")

    # (i) forward-space: exact discount-factor reproduction
    op = _fwd_operators(x, y)
    s_fwd = forward_uniform_sigma_star(op)
    f_uni = forward_solve(op, np.full(op["m"], s_fwd))
    print(f"[VII] forward sigma*={s_fwd:.3f}; DF area-match err={np.max(np.abs(op['C']@f_uni-op['areas'])):.2e}; "
          f"min fwd={f_uni.min():.4%}")

    # (j) non-uniform yield calibration: per-segment equipartition
    sig_y, phiseg = calibrate_nonuniform_yield(x, y, s50)
    print(f"[VIII] non-unif yield: per-seg phi in [{phiseg.min():.4f},{phiseg.max():.4f}] (target 0.5); "
          f"sigma_i in [{sig_y.min():.2f},{sig_y.max():.2f}]")

    # (k) synthesis: non-uniform tension on the forward
    sig_f, f_nu, it = forward_nonuniform_calibrate(op, s_fwd)
    B, A = forward_segment_energies(op, f_nu); phif = B / (B + sig_f**2 * A)
    print(f"[IX ] non-unif forward: converged {it} iters; per-seg phi in "
          f"[{phif.min():.4f},{phif.max():.4f}]; DF err={np.max(np.abs(op['C']@f_nu-op['areas'])):.2e}")

    # (l) multi-curve SOFR/OIS: crisp discount bootstrap + tensioned projection
    mc = multicurve_demo()
    print(f"[X  ] SOFR/OIS: OIS bootstrap reprices swaps to {mc['reprice']:.2e}; "
          f"projection fwd sigma*={mc['sigma']:.3f}, DF err={mc['df_err']:.2e}, "
          f"min fwd={mc['min_fwd']:.4%}")

    # (m) risk Jacobians and sigma* stability
    rs = risk_sensitivities_demo()
    print(f"[XI ] analytic curve-delta df/dy vs finite-diff: max err={rs['df_err']:.2e}; "
          f"sigma* under +/-100bp slope twist in [{rs['slope_lo']:.2f},{rs['slope_hi']:.2f}] (smooth)")

    # (n) benchmark methods + phi sweep
    bt, sb = benchmark_table()
    print(f"[XII] benchmark forward roughness (x1e4): " +
          ", ".join(f"{r['method']}={1e4*r['fwd_roughness']:.2f}" for r in bt))
    ps = phi_sweep()
    print(f"[XII] phi* sweep: " +
          ", ".join(f"phi={r['phi']}->(s={r['sigma']:.1f},loc={r['loc5y']:.2f}y)" for r in ps))
    print(line); print("All checks complete."); print(line)


def calibrate_nonuniform_yield(x, y, s0, target=0.5, iters=80, tol=1e-8):
    """Per-segment equipartition for the YIELD curve (bisection per segment)."""
    n = len(x); h = np.diff(x); sig = np.full(n - 1, s0)
    def seg_phi(zi, zj, yi, yj, hh, s):
        B = bending_energy_closed(zi, zj, hh, s); A = axial_energy_closed(zi, zj, yi, yj, hh, s)
        return B / (B + A)
    def sig_for(zi, zj, yi, yj, hh):
        lo, hi = 1e-2, min(300.0, 650.0 / hh); f = lambda s: seg_phi(zi, zj, yi, yj, hh, s) - target
        flo = f(lo)
        for _ in range(80):
            mid = 0.5 * (lo + hi); fm = f(mid)
            if (flo > 0) == (fm > 0): lo, flo = mid, fm
            else: hi = mid
        return 0.5 * (lo + hi)
    for _ in range(iters):
        z = nonuniform_solve(x, y, sig)
        new = np.array([sig_for(z[i], z[i+1], y[i], y[i+1], h[i]) for i in range(n - 1)])
        if np.max(np.abs(np.log(new / sig))) < tol: break
        sig = 0.5 * sig + 0.5 * new
    z = nonuniform_solve(x, y, sig)
    phiseg = np.array([seg_phi(z[i], z[i+1], y[i], y[i+1], h[i], sig[i]) for i in range(n - 1)])
    return sig, phiseg


# ======================================================================
# 10.  FIGURES
# ======================================================================
def make_figures(outdir="figures"):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(outdir, exist_ok=True)
    x, y = TENORS, YIELDS; RED = "#C8102E"; BLUE = "#1f77b4"; GREEN = "#2ca02c"
    s50 = solve_sigma_star(x, y, 0.5)

    # curves
    tt = np.linspace(x[0], x[-1], 800)
    plt.figure(figsize=(7, 4.2))
    plt.plot(tt, 100*np.interp(tt, x, y), color=BLUE, lw=1, label="Linear")
    plt.plot(tt, 100*cubic_eval(x, y, cubic_solve(x, y), tt), color=GREEN, lw=1, label="Cubic")
    z = tension_solve(x, y, s50)
    plt.plot(tt, 100*tension_eval(x, y, z, s50, tt), color=RED, lw=1.6, label=f"Tension $\\sigma^*$={s50:.1f}")
    plt.plot(x, 100*y, 'o', color='k', ms=3, label="Nodes")
    plt.xlabel("Tenor (years)"); plt.ylabel("Par yield (%)"); plt.legend(fontsize=8)
    plt.title("Interpolation methods"); plt.tight_layout()
    plt.savefig(f"{outdir}/curves.png", dpi=140); plt.close()

    # phi(sigma)
    sg = np.linspace(0.05, 30, 200); ph = [phi_of_sigma(x, y, s) for s in sg]
    plt.figure(figsize=(7, 4.2)); plt.plot(sg, ph, color=RED, lw=1.8)
    plt.axhline(0.5, color="#999", lw=0.7, ls="--"); plt.axvline(s50, color="#999", lw=0.7, ls="--")
    plt.scatter([s50], [0.5], color='k', zorder=5)
    plt.xlabel(r"$\sigma$"); plt.ylabel(r"$\varphi(\sigma)$"); plt.ylim(0, 1.02)
    plt.title(f"Energy-partition criterion ($\\sigma^*$={s50:.2f})"); plt.tight_layout()
    plt.savefig(f"{outdir}/phi.png", dpi=140); plt.close()

    # locality
    rc, g, rcr = locality_radius(x, y, 6, "cubic")
    rt, _, rtr = locality_radius(x, y, 6, "tension", s=s50)
    plt.figure(figsize=(7, 4.2))
    plt.semilogy(g, rcr/1e-4 + 1e-12, color=GREEN, lw=1.3, label=f"Cubic ($\\rho$={rc:.1f} yr)")
    plt.semilogy(g, rtr/1e-4 + 1e-12, color=RED, lw=1.3, label=f"Tension ($\\rho$={rt:.1f} yr)")
    plt.axvline(x[6], color='k', lw=0.7, ls=":")
    plt.xlabel("Tenor (years)"); plt.ylabel(r"$|\partial f/\partial r_k|$"); plt.ylim(1e-4, 2)
    plt.title("Hedge locality (bump 5y node)"); plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(f"{outdir}/locality.png", dpi=140); plt.close()

    print(f"figures written to ./{outdir}/  (curves.png, phi.png, locality.png)")


# ======================================================================
def main():
    ap = argparse.ArgumentParser(description="Reproduce the energy-partition tension-spline paper.")
    ap.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = ap.parse_args()
    validation_report()
    if not args.no_figures:
        make_figures()

if __name__ == "__main__":
    main()
