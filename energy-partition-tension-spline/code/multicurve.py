"""
Multi-curve SOFR/OIS demonstration (Sec. X of the paper).

Discount curve  : SOFR OIS, bootstrapped from par swap rates, held CRISP (exact).
Projection curve: a tenor-basis forward curve, TENSIONED by the energy-partition criterion.

SOFR par-swap levels are representative of the 22-May-2026 environment: anchored to
5Y SOFR ~ 4.25% (BlueGamma) and the Fed-funds target 3.50-3.75% (Mar-2026 FOMC), with a
negative long-end swap-spread shape (humped 20-30y). The tenor basis is illustrative.
Levels are representative; the construction is exact.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from tension_curve import (cubic_solve, cubic_eval, tension_solve, tension_eval,
                           solve_sigma_star, _fwd_operators, forward_solve,
                           forward_uniform_sigma_star, forward_segment_energies)

RED="#C8102E"; BLUE="#1f77b4"; GREEN="#2ca02c"; GREY="#888"

# ---- representative SOFR OIS par swap rates (annual, decimals) ----
sw_t = np.array([1.,2.,3.,5.,7.,10.,20.,30.])
sw_s = np.array([3.90,4.10,4.18,4.25,4.30,4.36,4.55,4.45])/100.0

# ---- 1. interpolate par rates to an annual grid, then bootstrap OIS discount factors ----
yrs = np.arange(1,31)
S   = np.interp(yrs, sw_t, sw_s)              # annual par swap rates
P   = np.zeros(len(yrs)); A = 0.0             # P[k]=DF(yrs[k]); A=running annuity
for k,Sn in enumerate(S):
    P[k] = (1.0 - Sn*A)/(1.0 + Sn)            # OIS par-swap bootstrap, tau=1
    A   += P[k]
# verify the bootstrap reprices every input swap exactly
annu = np.cumsum(P)
par_repriced = (1.0 - P)/annu
boot_err = np.max(np.abs(par_repriced - S))
print(f"OIS bootstrap reprices par swaps to max error {boot_err:.2e}")

z_ois = -np.log(P)/yrs                         # OIS continuously-compounded zero rates
# OIS instantaneous forward at annual nodes (backward diff of -lnP)
f_ois = np.empty(len(yrs)); f_ois[0] = -np.log(P[0])/1.0
f_ois[1:] = -(np.log(P[1:])-np.log(P[:-1]))/1.0

# ---- 2. projection curve: OIS zero + illustrative tenor basis (5->18 bp) ----
basis_t = (5 + 13*(sw_t-1)/29)/1e4            # bp at the swap tenors; illustrative
z_ois_t = -np.log(np.interp(sw_t, yrs, P))/sw_t   # OIS zero at swap tenors
z_proj_t = z_ois_t + basis_t

# ---- 3. tension the PROJECTION INSTANTANEOUS FORWARD directly (Sec. VII/IX) ----
#        area constraints reproduce the projection discount factors exactly.
op = _fwd_operators(sw_t, z_proj_t)
sig_fwd = forward_uniform_sigma_star(op, 0.5)
f_proj = forward_solve(op, np.full(op["m"], sig_fwd))
df_err = np.max(np.abs(op["C"]@f_proj - op["areas"]))
print(f"projection forward sigma*={sig_fwd:.3f}; projection DF area-match err={df_err:.2e}")
print(f"min projection fwd (tensioned) = {f_proj.min():.4%}  (>0 = no-arb)")
g = op["g"]

# cubic-in-yield-space projection forward, for contrast (the wiggly one)
zc = cubic_solve(sw_t, z_proj_t)
proj_cubic = cubic_eval(sw_t, z_proj_t, zc, g)
fwd_cubic  = proj_cubic + g*cubic_eval(sw_t, z_proj_t, zc, g, 1)
print(f"min projection fwd (cubic)     = {fwd_cubic.min():.4%}")

# OIS discount curve held crisp (exact DFs at all annual nodes by bootstrap)
opd = _fwd_operators(sw_t, z_ois_t)
f_ois_curve = forward_solve(opd, np.full(opd["m"], sig_fwd))

# ---- 4. figure ----
fig,(axL,axR)=plt.subplots(1,2,figsize=(10.6,4.2))
ttz=np.linspace(1,30,800)
axL.plot(yrs,100*z_ois,'-o',color=BLUE,lw=1.4,ms=2.5,label="OIS discount zero (crisp)")
axL.plot(sw_t,100*z_proj_t,'s',color=RED,ms=5,label="Projection zero nodes")
axL.set_xlabel("Tenor (years)"); axL.set_ylabel("Zero rate (%)")
axL.set_title("Two curves: OIS discount vs projection"); axL.legend(fontsize=8)

axR.plot(g,100*fwd_cubic,color=GREEN,lw=1.1,ls=":",label="Projection fwd (cubic-in-yield)")
axR.plot(g,100*f_proj,color=RED,lw=1.8,label=f"Projection fwd (tensioned $\\sigma^*$={sig_fwd:.1f})")
axR.plot(g,100*f_ois_curve,color=BLUE,lw=1.0,ls="--",label="OIS fwd (reference)")
axR.set_xlabel("Tenor (years)"); axR.set_ylabel("Instantaneous forward (%)")
axR.set_title("Projection forward: tensioned (smooth) vs cubic"); axR.legend(fontsize=8)
plt.tight_layout(); plt.savefig("fig_multicurve.png",dpi=150)
print("figure written: fig_multicurve.png")
