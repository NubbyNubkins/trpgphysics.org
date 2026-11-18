# Are the neutrinos and their hidden partners ladders like the others?
import numpy as np
nu=np.array([0.36e-3,(0.36**2+74.2)**.5*1e-3,(0.36**2+2510)**.5*1e-3])   # eV
le=np.array([0.51099895e6,105.6584e6,1776.93e6])                          # eV
N=le**2/nu
def shape(m):
    s=np.sqrt(np.array(m,float)); n=np.ones(3)/3**.5; par=s@n; perp=s-par*n
    op=np.degrees(np.arctan(np.linalg.norm(perp)/par)); b=np.array([0,0,1.])-n/3**.5; b/=np.linalg.norm(b)
    return op,np.degrees(np.arccos(perp@b/np.linalg.norm(perp)))
for k,m in (("neutrinos",nu),("leptons",le),("hidden partners",N)):
    op,tw=shape(m); print(f"{k:16s} steps x{m[1]/m[0]:8.1f} x{m[2]/m[1]:6.1f}  closing /{m[2]/m[0]:10.0f}  opening {op:5.2f}  twist {tw:5.2f}  twist/(54.74-op) {tw/(54.7356-op):.2f}")
print("hidden partners (GeV):",N/1e9)
print("check: hidden step = lepton step^2 / neutrino step:",(le[1]/le[0])**2/(nu[1]/nu[0]),(le[2]/le[1])**2/(nu[2]/nu[1]))
print("MMG 2.41e8 GeV vs hidden partners: between N1 and N2; sqrt(N1*N3) =",np.sqrt(N[0]*N[2])/1e9,"GeV")
print("\nDependence on the lightest neutrino (m1, meV): hidden ladder opening, twist, ratio; sqrt(N1 N3) GeV")
for m1 in (0.1,0.25,0.36,0.5,1.0,2.0):
    nu=np.array([m1,(m1**2+74.2)**.5,(m1**2+2510)**.5])*1e-3; N=le**2/nu; op,tw=shape(N); o2,t2=shape(nu)
    print(f"  m1 {m1:4.2f}: hidden {op:6.3f} {tw:6.3f} ratio {tw/(54.7356-op):.3f} | neutrino ratio {t2/(54.7356-o2):.3f} | sqrt(N1N3) {np.sqrt(N[0]*N[2])/1e9:.3e}")
op,tw=shape(le); print("lepton ratio",tw/(54.7356-op))
