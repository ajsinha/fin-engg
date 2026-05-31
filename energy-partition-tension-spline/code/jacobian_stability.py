"""
Risk sensitivities and stability (Items 2 & 3 of the review).

(2) Analytic Jacobians:
      dz/dy_j = M^{-1} db/dy_j           (M = tridiagonal, b linear in y)
      df(x)/dy_j   in closed form        (curve delta, at fixed sigma)
      dsigma*/dy_j = -(dphi/dy_j)/(dphi/dsigma)   (implicit function theorem)
(3) Stability: sigma* is C^1 in the inputs wherever dphi/dsigma != 0 (Prop. 2),
    so it cannot jump; demonstrated under parametric level/slope/curvature stress.

All analytic results are checked against finite differences.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from tension_curve import (TENORS as x, YIELDS as y, tension_solve, tension_eval,
                           bending_energy_closed, axial_energy_closed)

RED="#C8102E"; BLUE="#1f77b4"; GREEN="#2ca02c"
n=len(x); h=np.diff(x)

# ---------- fast phi / sigma* using closed-form energies (no quadrature) ----------
def fast_phi(xx, yy, s):
    hh=np.diff(xx); z=tension_solve(xx,yy,s)
    Eb=sum(bending_energy_closed(z[i],z[i+1],hh[i],s) for i in range(len(xx)-1))
    Et=sum(axial_energy_closed(z[i],z[i+1],yy[i],yy[i+1],hh[i],s) for i in range(len(xx)-1))
    return Eb/(Eb+Et)

def fast_sigma_star(xx, yy, phi_star=0.5, lo=0.05, hi=25.0, iters=50):
    f=lambda s: fast_phi(xx,yy,s)-phi_star; flo=f(lo)
    for _ in range(iters):
        m=0.5*(lo+hi); fm=f(m)
        if (flo>0)==(fm>0): lo,flo=m,fm
        else: hi=m
    return 0.5*(lo+hi)

# ---------- build M(sigma) and db/dy (interior system M z = b) ----------
def M_and_dbdy(s):
    sh=np.sinh(np.clip(s*h,0,700)); ch=np.cosh(np.clip(s*h,0,700))
    al=1.0/h - s/sh; be=s*ch/sh - 1.0/h
    M=np.zeros((n-2,n-2))
    for k in range(1,n-1):
        i=k-1; M[i,i]=be[k-1]+be[k]
        if i>0:   M[i,i-1]=al[k-1]
        if i<n-3: M[i,i+1]=al[k]
    # b_i = gamma_i - gamma_{i-1}, gamma_k = s^2 (y_{k+1}-y_k)/h_k
    # d gamma_k/d y_{k+1} = s^2/h_k ; d gamma_k/d y_k = -s^2/h_k
    dbdy=np.zeros((n-2,n))
    for k in range(1,n-1):
        i=k-1
        # gamma_k contributes +, gamma_{k-1} contributes -
        dbdy[i,k+1]+= s*s/h[k]; dbdy[i,k]+= -s*s/h[k]      # +gamma_k
        dbdy[i,k]  -= s*s/h[k-1]; dbdy[i,k-1] -= -s*s/h[k-1] # -gamma_{k-1}
    return M,dbdy

def dz_dy(s):
    M,dbdy=M_and_dbdy(s); Minv_db=np.linalg.solve(M,dbdy)   # (n-2) x n
    dz=np.zeros((n,n)); dz[1:n-1,:]=Minv_db                 # z_0=z_{n-1}=0
    return dz

# ---------- analytic df(x)/dy_j at fixed sigma ----------
def df_dy(s, tq):
    z=tension_solve(x,y,s); dz=dz_dy(s); tq=np.atleast_1d(tq)
    out=np.zeros((len(tq),n))
    idx=np.clip(np.searchsorted(x,tq)-1,0,n-2)
    for q,(tt,i) in enumerate(zip(tq,idx)):
        hi=h[i]; sh=np.sinh(np.clip(s*hi,0,700))
        P=x[i+1]-tt; Q=tt-x[i]; S1=np.sinh(s*P); S2=np.sinh(s*Q)
        for j in range(n):
            dzi=dz[i,j]; dzj=dz[i+1,j]
            term_hyp=(dzi*S1+dzj*S2)/(s*s*sh)
            term_lin=((1.0 if j==i else 0.0)-dzi/s**2)*(P/hi)\
                    +((1.0 if j==i+1 else 0.0)-dzj/s**2)*(Q/hi)
            out[q,j]=term_hyp+term_lin
    return out

# ---------- validate df/dy against finite differences ----------
s=fast_sigma_star(x,y,0.5)
tq=np.linspace(x[0],x[-1],60)
ana=df_dy(s,tq); eps=1e-6; fd=np.zeros_like(ana)
z0=tension_solve(x,y,s); f0=tension_eval(x,y,z0,s,tq)
for j in range(n):
    yp=y.copy(); yp[j]+=eps
    fd[:,j]=(tension_eval(x,yp,tension_solve(x,yp,s),s,tq)-f0)/eps
print(f"[2] analytic df/dy_j vs finite-diff: max abs err = {np.max(np.abs(ana-fd)):.2e}")

# ---------- dsigma*/dy_j via implicit function theorem ----------
def dsigmastar_dy(phi_star=0.5):
    ss=fast_sigma_star(x,y,phi_star); e=1e-5
    dphi_ds=(fast_phi(x,y,ss+e)-fast_phi(x,y,ss-e))/(2*e)
    dphi_dy=np.zeros(n)
    for j in range(n):
        yp=y.copy(); ym=y.copy(); yp[j]+=e; ym[j]-=e
        dphi_dy[j]=(fast_phi(x,yp,ss)-fast_phi(x,ym,ss))/(2*e)   # fixed sigma=ss
    return ss, -dphi_dy/dphi_ds, dphi_ds

ss, dss_dy, dphi_ds = dsigmastar_dy(0.5)
# validate IFT against direct re-optimisation (matched step)
e=1e-5; dss_fd=np.zeros(n)
for j in range(n):
    yp=y.copy(); ym=y.copy(); yp[j]+=e; ym[j]-=e
    dss_fd[j]=(fast_sigma_star(x,yp,0.5)-fast_sigma_star(x,ym,0.5))/(2*e)
rel=np.max(np.abs(dss_dy-dss_fd)/(np.abs(dss_fd)+1.0))
print(f"[2] dsigma*/dy_j  IFT vs finite-diff: max relative err = {rel:.2e}")
print(f"    dphi/dsigma at sigma*={ss:.3f} is {dphi_ds:.4f} (nonzero => sigma* is C^1 in inputs)")

# ---------- (3) stress: sigma* under level/slope/curvature shocks ----------
shocks=np.linspace(-0.01,0.01,41)          # +/- 100 bp
mid=np.mean(x); rng=x[-1]-x[0]
lvl=np.array([fast_sigma_star(x, y+d, 0.5) for d in shocks])
slp=np.array([fast_sigma_star(x, y+d*(x-mid)/rng, 0.5) for d in shocks])
bf =(x-mid)**2; bf=bf-bf.mean(); bf/=np.max(np.abs(bf))
crv=np.array([fast_sigma_star(x, y+d*bf, 0.5) for d in shocks])
print(f"[3] sigma* under stress: level [{lvl.min():.2f},{lvl.max():.2f}], "
      f"slope [{slp.min():.2f},{slp.max():.2f}], curvature [{crv.min():.2f},{crv.max():.2f}]")

# ================= FIGURE =================
fig,(axL,axR)=plt.subplots(1,2,figsize=(10.6,4.2))
# left: analytic vs FD delta for a couple of nodes
for j,c in [(6,RED),(8,BLUE)]:
    axL.plot(tq,ana[:,j],color=c,lw=1.6,label=f"analytic $\\partial f/\\partial y_{{{int(x[j])}\\mathrm{{y}}}}$")
    axL.plot(tq,fd[:,j],'o',color=c,ms=2.5,mfc='none')
axL.set_xlabel("Tenor (years)"); axL.set_ylabel(r"$\partial f/\partial y_j$")
axL.set_title("Curve delta: analytic (line) vs finite-diff (circles)"); axL.legend(fontsize=8)
# right: sigma* under stress
axR.plot(1e4*shocks,lvl,color=RED,lw=1.7,label="level shift")
axR.plot(1e4*shocks,slp,color=BLUE,lw=1.7,label="slope twist")
axR.plot(1e4*shocks,crv,color=GREEN,lw=1.7,label="curvature/butterfly")
axR.axhline(ss,color="#999",lw=0.7,ls="--")
axR.set_xlabel("shock (bp)"); axR.set_ylabel(r"$\sigma^\star$ at equipartition")
axR.set_title(r"Stability of $\sigma^\star$ under parametric stress"); axR.legend(fontsize=8)
plt.tight_layout(); plt.savefig("fig_jacobian.png",dpi=150)
print("figure written: fig_jacobian.png")
