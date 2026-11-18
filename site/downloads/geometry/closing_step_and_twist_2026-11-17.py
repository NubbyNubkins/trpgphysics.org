# Close each ladder (top rung back to bottom rung); are the uneven steps the twist? is the twist a steady torsion?
import numpy as np
F={"up":[2.16e-3,1.273,172.57],"down":[4.70e-3,0.0935,4.183],"leptons":[0.51099895e-3,0.1056584,1.77693],
   "neutrinos (m1 .36 meV)":[0.36,(0.36**2+74.2)**.5,(0.36**2+2510)**.5]}
def shape(m):
    s=np.sqrt(np.array(m,float)); n=np.ones(3)/3**.5; par=s@n; perp=s-par*n
    op=np.degrees(np.arctan(np.linalg.norm(perp)/par)); b=np.array([0,0,1.])-n/3**.5; b/=np.linalg.norm(b)
    tw=np.degrees(np.arccos(perp@b/np.linalg.norm(perp))); return op,tw
def masses(op,tw):
    e=2**.5*np.tan(np.radians(op)); s=np.array([1+e*np.cos(np.radians(tw+120*k)) for k in (1,2,0)]); m=np.sort(s**2); return m
print(f"{'ladder':24s} step1    step2    closing step      opening twist  twist/(54.74-opening)  steps if twist were 0")
for k,m in F.items():
    a,b=m[1]/m[0],m[2]/m[1]; op,tw=shape(m); m0=masses(op,0.0)
    print(f"{k:24s} x{a:7.1f} x{b:7.1f}  divide by {a*b:9.0f}  {op:6.2f} {tw:6.2f}  {tw/(54.7356-op):.3f}              x{m0[1]/m0[0]:.2f} then x{m0[2]/m0[1]:.1f}")
print("\nlog steps add to zero round the triangle: e.g. up +2.770 +2.132 -4.902 = 0")
print("successive twists (nu, lepton, down, up):",[round(shape(m)[1],2) for m in (F['neutrinos (m1 .36 meV)'],F['leptons'],F['down'],F['up'])])
