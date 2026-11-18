# Moebius flip and the cone repeated outward
import numpy as np
le=np.array([0.51099895e6,105.6584e6,1776.93e6]); up=np.array([2.16e6,1273e6,172570e6]); dn=np.array([4.70e6,93.5e6,4183e6])
def nus(m1): return np.array([m1,(m1**2+74.2)**.5,(m1**2+2510)**.5])*1e-3
GRAIN=8.63e27
def shape(m):
    s=np.sqrt(m); A=s.mean(); c=s/A-1; z=(2/3)*sum(c[i]*np.exp(-1j*np.radians(120*i)) for i in range(3))
    e=abs(z); th=np.degrees(np.angle(z)); return e,th+240-360,np.degrees(np.arctan(e/np.sqrt(2)))
nu=nus(0.36); N=le**2/nu; X4=N**2/le; X0=nu**2/le
print("(1) Sign flips of each ladder's wave in one full turn (Moebius needs an ODD number)")
for n,m in (("neutrinos",nu),("leptons",le),("hidden",N),("down",dn),("up",up)):
    e,tw,op=shape(m); print(f"  {n:10s} wave height/average e={e:.3f} -> {'2 flips (over and back)' if e>1 else '0 flips'}; stretch on the other face = {2*np.degrees(np.arccos(1/e)):.1f} deg of the 360")
print("  spin 1/2: one full turn multiplies the wave by -1 (1 flip per turn = Moebius); two turns restore it. Measured: Rauch 1975, Werner 1975 (neutrons).")
print("\n(2) The cone repeated: each copy = mirror of the one before through the one between (per rung)")
for n,m in (("copy 0 (inward)",X0),("1 neutrinos",nu),("2 leptons",le),("3 hidden",N),("4 next out",X4)):
    e,tw,op=shape(m); gm=np.exp(np.log(m).mean())
    print(f"  {n:16s} rungs {m[0]:.3g} {m[1]:.3g} {m[2]:.3g} eV | middle {gm:.3g} eV | steps x{m[1]/m[0]:.4g} x{m[2]/m[1]:.4g} | twist {tw:6.2f} opening {op:5.2f} gap-to-edge {54.7356-op:5.2f}")
gms=[np.exp(np.log(m).mean()) for m in (nu,le,N,X4)]; print("  equal jump between copies: x%.3g"%(gms[1]/gms[0]))
print("  grain = %.3g eV; copy 4 rungs / grain ="%GRAIN, X4/GRAIN)
print("  copy 4 for other lightest neutrinos (rung 1 / grain):",[(m1,float((le[0]**2/nus(m1)[0])**2/le[0]/GRAIN)) for m1 in (0.1,0.36,1.0)])
tw=[shape(m)[1] for m in (nu,le,N,X4)]; gp=[54.7356-shape(m)[2] for m in (nu,le,N,X4)]
print("  twist shrink per copy:",np.round(np.array(tw[1:])/np.array(tw[:-1]),3),"  gap-to-edge shrink:",np.round(np.array(gp[1:])/np.array(gp[:-1]),3))
print("  compare up quarks: twist %.2f gap %.2f ; down: twist %.2f gap %.2f"%(shape(up)[1],54.7356-shape(up)[2],shape(dn)[1],54.7356-shape(dn)[2]))
