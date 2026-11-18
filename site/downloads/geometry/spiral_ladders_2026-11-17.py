# Ladders as spiral waves, not triangles (rotation matrix + shrink each step)
import numpy as np
le=np.array([0.51099895e6,105.6584e6,1776.93e6]); up=np.array([2.16e6,1273e6,172570e6]); dn=np.array([4.70e6,93.5e6,4183e6])
def nus(m1): return np.array([m1,(m1**2+74.2)**.5,(m1**2+2510)**.5])*1e-3
nu=nus(0.36); N=le**2/nu
L={"neutrinos":nu,"leptons":le,"hidden":N,"down":dn,"up":up}
print("(a) Each ladder as ONE triangle turned by the rotation formula and scaled")
def cone(m):
    s=np.sqrt(m); A=s.mean(); c=s/A-1            # c_i = e cos(delta + 120 i), i = heaviest..: fit
    # complex amplitude
    z=(2/3)*sum(c[i]*np.exp(-1j*np.radians(120*i)) for i in range(3)); return A,abs(z),np.degrees(np.angle(z))
for k,m in L.items():
    A,e,d=cone(m[::-1]); tri=np.array([[np.cos(np.radians(120*i)),np.sin(np.radians(120*i))] for i in range(3)])
    th=np.radians(d); R=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]]); rot=(R@tri.T).T   # x'=x cos - y sin ; y'=x sin + y cos
    back=(A*(1+e*rot[:,0]))**2
    print(f"  {k:10s} size e={e:.4f}  turn={abs(d):6.2f} deg  rebuilt/actual = {np.round(back/m[::-1],6)}")
print("\n(b) Is one ladder a pure spiral (same turn, same factor each step)?  shrink k = ln(step2)/ln(step1); pure log spiral needs k=1")
def k(m): return np.log(m[2]/m[1])/np.log(m[1]/m[0])
for n,m in L.items(): print(f"  {n:10s} steps x{m[1]/m[0]:9.1f} then x{m[2]/m[1]:6.1f}   k = {k(m):.3f}")
upZ=np.array([1.27e6,0.619e9,171.7e9]); dnZ=np.array([2.90e6,55e6,2.89e9])
print(f"  at the Z footing: up k = {k(upZ):.3f}, down k = {k(dnZ):.3f}")
print("\n(c) lightest neutrino for which the neutrino spiral shrinks exactly like the lepton spiral")
kl=k(le)
for m1 in (0.0,0.1,0.25,0.36,0.5,0.75,1.0,2.0,5.0): print(f"  m1 {m1:4.2f} meV: neutrino k {k(nus(max(m1,1e-6))) if m1>0 else float('nan'):.3f}  hidden k {k(le**2/nus(max(m1,1e-6))):.3f}   (lepton {kl:.3f})")
lo,hi=0.01,5.0
for _ in range(80):
    mid=(lo+hi)/2
    if k(nus(mid))<kl: lo=mid
    else: hi=mid
print(f"  equal at m1 = {mid:.3f} meV")
# chance: m1 log-uniform 0.01..30 meV, how often |k_nu - k_lep| < observed at 0.36
obs=abs(k(nus(0.36))-kl); g=np.exp(np.random.default_rng(1).uniform(np.log(0.01),np.log(30),200000)); kk=np.array([k(nus(x)) for x in g[:20000]])
print(f"  chance a random m1 (0.01-30 meV) lands this close: {np.mean(abs(kk-kl)<obs)*100:.0f}%")
print("\n(d) Across ladders: one spiral through all of them? (turn angle vs size, view from above)")
rows=[]
for n,m in L.items():
    A,e,d=cone(m[::-1]); op=np.degrees(np.arctan(e/np.sqrt(2))); rows.append((n,abs(d),e,op))
rows.sort(key=lambda r:r[1])
for n,d,e,op in rows: print(f"  {n:10s} turn {d:6.2f}  size e {e:.3f}  opening {op:5.2f}  Archimedean: turn/(54.74-opening) = {d/(54.7356-op):.2f}   log: ln(2/e)/turn = {np.log(2/e)/d:.4f}")
t=[r[1] for r in rows]; print("  turn ratios ladder to ladder:",np.round(np.array(t[1:])/np.array(t[:-1]),2))
