#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
locality_targeting.py -- PROTOTYPE / proof-of-concept.

Question: can a desk specify a per-node hedge-locality PROFILE (e.g. tight at 5y,
looser at 10y) and have the model return the per-segment tension vector that best
meets it while keeping the forward curve smooth?

This is a genuine *optimization*, not a clean inversion, because the spline is C^2
and globally coupled: the locality at node k depends on the WHOLE tension vector,
not just the adjacent segment. So we minimize

    J(sigma) = sum_k w_k ( rho_k(sigma) - rho_k^target )^2   +   lam * R(sigma)

over the per-segment tension vector sigma >= 0, where rho_k is the hedge-locality
radius at node k and R is the forward roughness (smoothness budget, weight lam).

What this script reports, ON THE REAL 22-May-2026 curve (nothing simulated):
  0. a per-node FEASIBILITY FLOOR: the tightest locality each node can reach, set
     by knot spacing -- tension compresses locality only down to ~ the local mesh,
     so targets below the floor are unreachable (an honest, computable limit),
  1. a gradient sanity check: two finite-difference step sizes agree,
  2. an optimization run toward a desk-specified locality profile, reporting the
     achieved-vs-requested locality (FEASIBILITY reported, not assumed),
  3. a smoothness/locality FRONTIER as the budget weight lam varies,
  4. a convexity probe of R(sigma) along random directions (honesty check).

KEY FINDING (real data): per-node locality targeting works well in the DENSE region
of the curve (front/belly, where desks hold the most instruments), and is MESH-LIMITED
at the sparse long end (20-30y), where no tension delivers tight locality. The model
reports the floor rather than over-promising.

Honest caveats are printed inline. Requires numpy, scipy (optional), matplotlib,
and tension_curve.py.

Usage:
    python locality_targeting.py
    python locality_targeting.py --targets 5:0.4 10:1.2 20:2.0
"""
from __future__ import annotations
import argparse
import numpy as np

import tension_curve as tc

X = tc.TENORS.astype(float)
Y = tc.YIELDS.astype(float)
M = len(X) - 1                         # number of segments
GRID = np.linspace(X[0], X[-1], 600)
LGRID = np.linspace(X[0], X[-1], 300)  # coarser grid for locality ratio integrals


# ----------------------------------------------------------------------
# curve + metrics under a per-segment tension VECTOR
# ----------------------------------------------------------------------
def curve_zero(sig):
    z = tc.nonuniform_solve(X, Y, sig)
    return tc.nonuniform_eval(X, Y, z, sig, GRID)

def fwd_roughness(sig):
    zg = curve_zero(sig)
    f = np.gradient(GRID * zg, GRID)
    return float(np.trapezoid(np.gradient(f, GRID) ** 2, GRID))

def locality_at(sig, k, eps=1e-4):
    """Hedge-locality radius at node k under tension vector sig (bumped-node response)."""
    z0 = tc.nonuniform_solve(X, Y, sig)
    f0 = tc.nonuniform_eval(X, Y, z0, sig, LGRID)
    yp = Y.copy(); yp[k] += eps
    z1 = tc.nonuniform_solve(X, yp, sig)
    f1 = tc.nonuniform_eval(X, yp, z1, sig, LGRID)
    resp = np.abs(f1 - f0)
    denom = np.trapezoid(resp, LGRID)
    if denom <= 0:
        return np.nan
    return float(np.trapezoid(np.abs(LGRID - X[k]) * resp, LGRID) / denom)

def locality_vector(sig, nodes):
    return np.array([locality_at(sig, k) for k in nodes])


# ----------------------------------------------------------------------
# objective and its finite-difference gradient
# (locality is a coupled functional; we use finite-difference gradients here,
#  which is honest for a prototype. Analytic df/dy exists in tension_curve;
#  extending it to drho/dsigma is the production follow-up.)
# ----------------------------------------------------------------------
def objective(sig, nodes, targets, weights, lam):
    rho = locality_vector(sig, nodes)
    fit = np.sum(weights * (rho - targets) ** 2)
    return fit + lam * fwd_roughness(sig), rho

def grad_fd(fn, sig, h=1e-4):
    g = np.zeros_like(sig); f0 = fn(sig)
    for j in range(len(sig)):
        sp = sig.copy(); sp[j] += h
        g[j] = (fn(sp) - f0) / h
    return g


# ----------------------------------------------------------------------
# optimizer: projected gradient descent on log(sigma) (keeps sigma > 0)
# ----------------------------------------------------------------------
def optimize(nodes, targets, weights, lam, s0=2.0, iters=50, lr=0.15, tol=1e-6, gtol=5e-3, verbose=False):
    u = np.log(np.full(M, s0))                      # u = log sigma, unconstrained
    Jfn = lambda uu: objective(np.exp(uu), nodes, targets, weights, lam)[0]
    best = (np.inf, u.copy()); prev = np.inf
    for it in range(iters):
        g = grad_fd(Jfn, u, h=1e-3)
        gn = np.linalg.norm(g)
        if gn > 0:
            u = u - lr * g / max(gn, 1e-8)          # normalized step (robust to scale)
        J = Jfn(u)
        if J < best[0]:
            best = (J, u.copy())
        if abs(prev - J) < tol or gn < gtol:        # converged (value flat or small gradient)
            break
        prev = J
        if verbose and it % 25 == 0:
            print(f"   iter {it:3d}  J={J:.6e}  |grad|={gn:.2e}")
    sig = np.exp(best[1])
    return sig, best[0]


# ----------------------------------------------------------------------
def parse_targets(items):
    """['5:0.4','10:1.2'] -> (node_indices, target_radii)."""
    nodes, tgt = [], []
    for it in items:
        ten, rho = it.split(":")
        k = int(np.argmin(np.abs(X - float(ten))))
        nodes.append(k); tgt.append(float(rho))
    return np.array(nodes), np.array(tgt, float)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", nargs="+", default=["5:0.9", "7:1.3", "10:2.0"],
                    help="node:radius pairs, e.g. 5:0.9 7:1.3 10:2.0")
    ap.add_argument("--lam", type=float, default=2.0, help="smoothness budget weight")
    ap.add_argument("--export", action="store_true",
                    help="write CSVs (floor, achieved-vs-target, lambda frontier) for charting/upload")
    args = ap.parse_args()

    nodes, targets = parse_targets(args.targets)
    weights = np.ones(len(nodes))
    names = [f"{X[k]:g}y" for k in nodes]
    print(f"curve: real 22-May-2026 U.S. Treasury ({M} segments)")
    print(f"desk locality targets: " + ", ".join(f"{n}->{t:.2f}y" for n, t in zip(names, targets)))

    # baseline: global equipartition tension (the paper's default)
    s_glob = tc.solve_sigma_star(X, Y, 0.5)
    base_rho = locality_vector(np.full(M, s_glob), nodes)
    print(f"\nbaseline (global sigma*={s_glob:.2f}) locality: " +
          ", ".join(f"{n}={r:.2f}y" for n, r in zip(names, base_rho)))
    print(f"baseline forward roughness x1e4 = {1e4*fwd_roughness(np.full(M, s_glob)):.2f}")

    # feasibility floor: tightest locality each node can reach (sigma -> large).
    # Tension compresses locality only down to ~ the local knot spacing; beyond that
    # the mesh, not the tension, controls it. Report this BEFORE optimizing.
    floor = locality_vector(np.full(M, 40.0), nodes)
    print("\n[0] feasibility floor (tightest achievable locality, sigma->large):")
    for n, fl, t in zip(names, floor, targets):
        ok = "reachable" if t >= fl - 0.05 else f"BELOW FLOOR ({fl:.2f}y) -> not achievable"
        print(f"      {n}: floor={fl:.2f}y, target={t:.2f}y  [{ok}]")
    print("    (locality cannot be driven below the knot-spacing floor by tension alone.)")
    export_rows = {"floor": [(names[i], float(floor[i]), float(targets[i])) for i in range(len(nodes))]}

    # 1) gradient sanity check at a random positive sigma
    print("\n[1] gradient check (objective dJ/d log-sigma, two FD step sizes):")
    rng = np.random.default_rng(0); u0 = np.log(rng.uniform(0.5, 4.0, M))
    Jfn = lambda uu: objective(np.exp(uu), nodes, targets, weights, args.lam)[0]
    g1 = grad_fd(Jfn, u0, 1e-3); g2 = grad_fd(Jfn, u0, 5e-4)
    rel = np.max(np.abs(g1 - g2) / (np.abs(g2) + 1e-9))
    print(f"    max rel. difference between FD steps = {rel:.2e}  (small => gradient is well-defined)")

    # 2) optimize to the target profile
    print("\n[2] optimizing tension vector to the locality profile...")
    sig, J = optimize(nodes, targets, weights, args.lam, verbose=True)
    ach = locality_vector(sig, nodes)
    print(f"\n    achieved locality: " + ", ".join(f"{n}={r:.2f}y (req {t:.2f})"
          for n, r, t in zip(names, ach, targets)))
    print(f"    residual |achieved-target| = {np.abs(ach-targets)}")
    print(f"    forward roughness x1e4 = {1e4*fwd_roughness(sig):.2f}  (baseline {1e4*fwd_roughness(np.full(M,s_glob)):.2f})")
    print(f"    per-segment tension sigma = {np.round(sig,2)}")
    feasible = np.all(np.abs(ach - targets) < 0.15)
    print(f"    => targets {'MET within 0.15y (feasible)' if feasible else 'NOT all met: best-achievable fit (some targets conflict / out of range)'}")
    export_rows["fit"] = [(names[i], float(targets[i]), float(floor[i]), float(ach[i])) for i in range(len(nodes))]

    # 3) frontier: sweep the smoothness weight lam
    print("\n[3] smoothness/locality frontier (sweeping budget weight lam):")
    print(f"    {'lam':>8s} {'rough x1e4':>11s} {'mean|rho-tgt|':>13s}")
    frontier = []
    for lam in [0.0, 4.0, 32.0]:
        s, _ = optimize(nodes, targets, weights, lam, iters=60)
        r = locality_vector(s, nodes)
        rough = 1e4 * fwd_roughness(s); fit = float(np.mean(np.abs(r - targets)))
        frontier.append((float(lam), rough, fit))
        print(f"    {lam:8.1f} {rough:11.2f} {fit:13.3f}")
    print("    (higher lam -> smoother forward but looser fit to the locality profile:")
    print("     this is the achievable trade-off the desk chooses on.)")

    # 4) convexity probe of R(sigma)
    print("\n[4] convexity probe of forward roughness R(sigma) along random directions:")
    s_mid = np.full(M, 2.0); bad = 0; ndir = 4
    for _ in range(ndir):
        d = rng.normal(size=M); d /= np.linalg.norm(d)
        ts = np.linspace(-0.6, 0.6, 9)
        vals = [fwd_roughness(np.exp(np.log(s_mid) + a*d)) for a in ts]
        sec = np.diff(vals, 2)
        if np.any(sec < -1e-6): bad += 1
    print(f"    non-convex directions: {bad}/{ndir}  "
          f"({'R looks locally convex -> single optimum likely' if bad==0 else 'R is non-convex -> use multistart; report best-achievable'})")

    if args.export:
        import csv
        with open("locality_targeting_results.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["section", "node", "target", "floor", "achieved"])
            for n, t, fl, a in export_rows["fit"]:
                w.writerow(["fit", n, t, fl, a])
            w.writerow([]); w.writerow(["section", "lambda", "roughness_x1e4", "mean_abs_resid"])
            for lam, rg, ft in frontier:
                w.writerow(["frontier", lam, rg, ft])
        print("\n[export] wrote locality_targeting_results.csv  (upload this for the chart/subsection)")

    print("\nHONEST READ: this is an operational optimization layer on top of the model,")
    print("not new clean theory. Locality is globally coupled, so targets may conflict;")
    print("the script reports best-achievable fit + the smoothness trade-off rather than")
    print("pretending every profile is reachable.")


if __name__ == "__main__":
    main()
