# What could make both Hubble values right without anything extra at the beginning?
# Test: atoms forming slightly earlier (ruler shorter). (1) clumpy matter when atoms formed (regions of different resistance); (2) a slightly stronger electron grip then.
# Simple three-level hydrogen atom (Peebles 1968), helium ignored. Used only for DIFFERENCES from its own baseline.
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
c=2.99792458e10; kB=8.617333e-5; hbar=6.582120e-16; mec2=510998.95; sigT0=6.6524587e-25; mH=1.6735e-24; Mpc=3.0857e24
T0=2.7255; wb=0.02237; wm=0.1430; wr=4.183e-5; Yp=0.245
nH0=(1-Yp)*wb*1.8788e-29/mH
def Hz(z,h): return 100*h*1e5/Mpc*np.sqrt(wr/h**2*(1+z)**4+wm/h**2*(1+z)**3+(1-(wm+wr)/h**2))
def make(d=0.0):
    s=1+d; E1,E2=13.6057*s,3.4014*s; E21=E1-E2; lam=1.215670e-5/s; Lam=8.2246*s; sig=sigT0/s**2
    def alpha(T): t=T/s/1e4; return 1.14*4.309e-13*t**-0.6166/(1+0.6703*t**0.53)/s**2
    def beta(T): return alpha(T)*(mec2*s*kB*T/(2*np.pi))**1.5/(hbar*c)**3*np.exp(-E2/(kB*T))
    return E1,E2,E21,lam,Lam,sig,alpha,beta
def xe_zone(Delta,h,d,zs):
    E1,E2,E21,lam,Lam,sig,alpha,beta=make(d)
    def saha(z):
        T=T0*(1+z); n=nH0*Delta*(1+z)**3; A=(mec2*(1+d)*kB*T/(2*np.pi))**1.5/(hbar*c)**3*np.exp(-E1/(kB*T))/n; return (-A+np.sqrt(A*A+4*A))/2
    def rhs(z,x):
        x=min(max(x[0],1e-10),1.0); T=T0*(1+z); n=nH0*Delta*(1+z)**3; H=Hz(z,h); K=lam**3/(8*np.pi*H); n1=(1-x)*n
        b=beta(T); C=(1+K*Lam*n1)/(1+K*(Lam+b)*n1)
        return [C*(alpha(T)*x*x*n-b*(1-x)*np.exp(-E21/(kB*T)))/(H*(1+z))]
    z0=1650*(1+d); sol=solve_ivp(rhs,[z0,zs[-1]],[saha(z0)],t_eval=zs[zs<=z0][::1],method="Radau",rtol=1e-6,atol=1e-10)
    x=np.ones_like(zs); x[zs<=z0]=sol.y[0]; x[zs>z0]=[min(1,saha(z)) for z in zs[zs>z0]]; return x
def zstar(zones,h,d=0.0):
    zs=np.linspace(2200,200,4001); E1,E2,E21,lam,Lam,sig,alpha,beta=make(d)
    xbar=sum(fv*D*xe_zone(D,h,d,zs) for fv,D in zones)      # mean electron density / mean n_H
    dt=xbar*nH0*(1+zs)**2*sig*c/Hz(zs,h)                   # dtau/dz
    tau=np.cumsum(dt[::-1])[::-1]*abs(zs[1]-zs[0]); tau=tau-tau[0]+0  # from z upward... recompute properly below
    tau=np.array([np.trapezoid(dt[i:][::-1],zs[i:][::-1]) for i in range(len(zs))])
    g=dt*np.exp(-tau); return zs[np.argmax(g)], zs[np.argmin(abs(tau-1))]
def rs(z,h): R=lambda zz: 3*wb/(4*2.4729e-5)/(1+zz); return quad(lambda zz: c/np.sqrt(3*(1+R(zz)))/Hz(zz,h),z,1e7,limit=400)[0]/Mpc
def DA(z,h): return quad(lambda zz: c/Hz(zz,h),0,z,limit=400)[0]/Mpc
base=[(1.0,1.0)]; h0=0.674; zp0,zt0=zstar(base,h0); th0=rs(zp0,h0)/DA(zp0,h0)
print(f"baseline (smooth matter, H0 67.4): light set free at stretch {zp0:.0f} (accurate codes: 1090); ruler {rs(zp0,h0):.1f} Mpc (accurate: 144.4 at this moment); angle {100*th0:.4f}")
def solve(zones,d=0.0,label=""):
    zp,_=zstar(zones,h0,d); r=rs(zp,h0)
    hh=brentq(lambda h: rs(zp,h)/DA(zp,h)-th0,0.55,0.95)
    b=sum(fv*D*D for fv,D in zones)-1
    print(f"  {label:58s} clumping b={b:.2f}: atoms form earlier by {100*(zp/zp0-1):+.2f}%; ruler {100*(r/rs(zp0,h0)-1):+.2f}%; same sky angle needs H0 = {100*hh:.1f}")
    return hh
print("\n(1) Clumpy matter when atoms formed (three equal-volume regions; mean density unchanged)")
solve(base,label="smooth")
for lo,hi,lab in ((0.6,1.4,"mild: densities 0.6 / 1 / 1.4"),(0.3,1.7,"medium: 0.3 / 1 / 1.7"),(0.1,1.9,"strong: 0.1 / 1 / 1.9 (Jedamzik & Pogosian 2020 model M1)"),(0.02,1.98,"extreme: 0.02 / 1 / 1.98")):
    solve([(1/3,lo),(1/3,1.0),(1/3,hi)],label=lab)
solve([(0.5,0.05),(0.25,1.0),(0.25,2.9)],label="voids and walls: half the volume nearly empty (0.05 / 1 / 2.9)")
print("\n(2) A slightly stronger electron grip then (heavier electron at atom formation, normal today)")
for d in (0.005,0.01,0.02,0.03): solve(base,d,label=f"electron heavier by {100*d:.1f}%")
print("\nneeded for H0 = 73.0: ruler about -7.7% if nothing else changes")
