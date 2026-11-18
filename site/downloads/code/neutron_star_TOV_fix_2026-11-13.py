import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
CONV=1.3234e-6; MSUN=1.4766
B=54.1; a=0.00042
def p_gem(e):
    if e<=4*B: return max(1e-4,0.12*e*(e/(4*B))**0.65)
    d=e-4*B; return d/3+a*d*d/(1+1.45*a*d)
egrid=np.logspace(-1,np.log10(5000),40000)
def make_inv(pf):
    pg=np.array([pf(e) for e in egrid]); pg=np.maximum.accumulate(pg+egrid*1e-12)
    return lambda p: np.interp(p,pg,egrid)
def tov(pf,einv,ec,psurf=1e-6):
    pc=pf(ec)
    def rhs(r,y):
        p,m=y
        if p<=psurf: return [0,0]
        e=einv(p)*CONV; pg=p*CONV
        return [-(e+pg)*(m+4*np.pi*r**3*pg)/(r*(r-2*m))/CONV, 4*np.pi*r*r*e]
    ev=lambda r,y: y[0]-psurf; ev.terminal=True
    r0=1e-4; s=solve_ivp(rhs,[r0,60],[pc,4/3*np.pi*r0**3*ec*CONV],events=ev,rtol=1e-8,atol=1e-12,max_step=0.02)
    return s.t[-1], s.y[1][-1]/MSUN
inv=make_inv(p_gem)
print("The earlier version's EOS, solved properly:")
rows=[]
for ec in [280,350,430,520,650,780,920,1080,1240,1380,1500,1650,2000,2500,3000]:
    R,M=tov(p_gem,inv,ec); rows.append((ec,M,R))
    cs2=(p_gem(ec+0.5)-p_gem(ec-0.5))
    print(f"eps_c={ec:5d}  M={M:.3f} Msun  R={R:.2f} km  cs2={cs2:.3f}")
# low-density part: p = 0.12 e (e/216)^0.65 at e=216 -> p=26 MeV/fm3 (nuclear saturation pressure ~ 0)
print("p at nuclear density e=150:",p_gem(150), " (real nuclear matter ~ 0-3 MeV/fm3)")


# ---- SLy reference check and bag-star scan ----
import numpy as np
# (uses tov(), p_gem, B, a defined above)
# low-density check for the earlier EOS
inv=make_inv(p_gem)
print("Earlier EOS, low central densities:")
for ec in [20,40,80,120,160,200,240]:
    R,M=tov(p_gem,inv,ec); print(f"  eps_c={ec} M={M:.3f} R={R:.1f}")
# SLy piecewise polytrope (Read et al. 2009), cgs, p in dyn/cm2 via K rho^G * c^2
c2=8.98755e20
crust=[(6.80110e-9,1.58425,2.44034e7),(1.06186e-6,1.28733,3.78358e11),(5.32697e1,0.62223,2.62780e12),(3.99874e-8,1.35692,None)]
G1,G2,G3=3.005,2.988,2.851; rho1,rho2=10**14.7,10**15.0; p1=10**34.384/c2
K1=p1/rho1**G1
K4c,G4c=crust[3][:2]; rho0=(K4c/K1)**(1/(G1-G4c))
pieces=[]; lo=0
for K,G,hi in crust[:3]: pieces.append([K,G,lo,hi]); lo=hi
pieces.append([K4c,G4c,lo,rho0]); pieces.append([K1,G1,rho0,rho1])
K2=K1*rho1**(G1-G2); pieces.append([K2,G2,rho1,rho2])
K3=K2*rho2**(G2-G3); pieces.append([K3,G3,rho2,1e17])
# energy-density constants a_i for continuity
a=0; prev=None
for pc in pieces:
    K,G,lo,hi=pc
    if prev is None: pc.append(0.0)
    else:
        Kp,Gp,_,_,ap=prev; e_lo=(1+ap)*lo+Kp/(Gp-1)*lo**Gp
        pc.append(e_lo/lo-1-K/(G-1)*lo**(G-1))
    prev=pc
MEVFM3=1.7827e12  # g/cm3 per MeV/fm3
rr=np.logspace(3,16.5,60000); P=[];E=[]
for r in rr:
    for K,G,lo,hi,ai in pieces:
        if lo<=r<hi: P.append(K*r**G); E.append((1+ai)*r+K/(G-1)*r**G); break
P=np.array(P)/MEVFM3; E=np.array(E)/MEVFM3
p_sly=lambda e: np.interp(e,E,P); inv_s=lambda p: np.interp(p,P,E)
print("\nReal nuclear EOS (SLy):  p at e=150:",p_sly(150))
best=[]
for ec in np.linspace(200,3000,57):
    R,M=tov(p_sly,inv_s,ec,psurf=1e-12); best.append((ec,M,R))
best=np.array(best); i=best[:,1].argmax()
print(f"SLy max mass {best[i,1]:.3f} Msun at R={best[i,2]:.2f} km (eps_c={best[i,0]:.0f})")
j=np.argmin(abs(best[:i,1]-1.4)); print(f"SLy 1.4 Msun star: R={best[j,2]:.2f} km (M={best[j,1]:.3f})")
import numpy as np

hc=197.327
for Bv in [54.1,57.0,60.0,70.0]:
  for av in [0.0,0.00042]:
    pf=lambda e,Bv=Bv,av=av: max((e-4*Bv)/3 + av*(e-4*Bv)**2/(1+1.45*av*(e-4*Bv)),0.0)
    einv=lambda p,pf=pf,Bv=Bv: brentq(lambda e: pf(e)-p, 4*Bv*0.999999, 1e5)
    out=[]
    for ec in np.linspace(4*Bv*1.01,3500,60):
        R,M=tov(pf,einv,ec,psurf=1e-7); out.append((ec,M,R))
    out=np.array(out); i=out[:,1].argmax(); j=np.argmin(abs(out[:i,1]-1.4))
    print(f"B={Bv:5.1f} (B^1/4={(Bv*hc**3)**.25:.1f} MeV) wall={av}: Mmax={out[i,1]:.2f} R_Mmax={out[i,2]:.2f}  R1.4={out[j,2]:.2f}")
