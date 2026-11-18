# K2: real tidal Love number k2 and Lambda (Hinderer 2008 method), self-bound surface correction
# Validation: SLy gives Lambda(1.4)=299 (published ~300).
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
CONV=1.3234e-6; MSUN=1.4766
def k2_of(C,y):
    t1=8*C**5/5*(1-2*C)**2*(2+2*C*(y-1)-y)
    t2=2*C*(6-3*y+3*C*(5*y-8))+4*C**3*(13-11*y+C*(3*y-2)+2*C*C*(1+y))+3*(1-2*C)**2*(2-y+2*C*(y-1))*np.log(1-2*C)
    return t1/t2
def star(pf,dpde,einv,ec,es_surf=0.0,psurf=1e-9):
    pc=pf(ec)
    def rhs(r,Y):
        p,m,y=Y
        e=einv(max(p,psurf))*CONV; P=p*CONV
        f=1-2*m/r; cs2=dpde(einv(max(p,psurf)))
        dp=-(e+P)*(m+4*np.pi*r**3*P)/(r*r*f)
        Q=4*np.pi/f*(5*e+9*P+(e+P)/cs2)-6/(r*r*f)-(2*(m+4*np.pi*r**3*P)/(r*r*f))**2
        F=(1-4*np.pi*r*r*(e-P))/f
        dy=-(y*y+y*F+r*r*Q)/r
        return [dp/CONV,4*np.pi*r*r*e,dy]
    ev=lambda r,Y:Y[0]-psurf; ev.terminal=True
    r0=1e-4
    s=solve_ivp(rhs,[r0,60],[pc,4/3*np.pi*r0**3*ec*CONV,2.0],events=ev,rtol=1e-9,atol=1e-14,max_step=0.01)
    R=s.t[-1]; M=s.y[1][-1]; y=s.y[2][-1]
    y-=4*np.pi*R**3*es_surf*CONV/M   # surface density jump correction (self-bound star)
    C=M/R; k2=k2_of(C,y); Lam=2/3*k2/C**5
    return M/MSUN,R,k2,Lam
def seq(name,B,av):
    pf=lambda e: max((e-4*B)/3+av*(e-4*B)**2/(1+1.45*av*(e-4*B)),0.0)
    h=1e-3
    dpde=lambda e: (pf(e+h)-pf(max(e-h,4*B)))/(e+h-max(e-h,4*B))
    einv=lambda p: brentq(lambda e: pf(e)-p,4*B,1e5) if p>0 else 4*B
    rows=[star(pf,dpde,einv,ec,es_surf=4*B) for ec in np.geomspace(4*B*1.01,2500,70)]
    rows=np.array(rows); i=rows[:,0].argmax(); r=rows[:i+1]
    L14=float(np.exp(np.interp(1.4,r[:,0],np.log(r[:,3])))); k14=np.interp(1.4,r[:,0],r[:,2])  # log(Lambda) interpolated: Lambda ~ C^-5
    print(f"{name}: Mmax={rows[i,0]:.2f}  k2(1.4)={k14:.3f}  Lambda(1.4)={L14:.0f}")
    return rows
# Main run: tidal test (GW170817 needs Lambda(1.4) below about 580; heaviest star at least 2.0 Msun)
if __name__ == "__main__":
    for B, av in [(54.1, 0.0), (62.0, 0.0), (62.0, 0.0002), (66.0, 0.0002), (70.0, 0.0002), (70.0, 0.0003)]:
        seq(f"B={B} wall={av}", B, av)
