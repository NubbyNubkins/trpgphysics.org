# The open items revisited with the spiral / wave / timeline tools
import numpy as np
rng=np.random.default_rng(7); ln=np.log; lg=np.log10
le=[0.51099895,105.6584,1776.93]
def nus(m1): return [m1,(m1**2+74.2)**.5,(m1**2+2510.)**.5]   # meV; NuFIT splittings 7.42e-5, 2.51e-3 eV^2
def f(m): return ln(m[1]/m[0])/ln(m[2]/m[0])       # where the middle rung sits between the ends (0 = bottom, 1 = top), on the scale of tens
print("(1) Where the middle rung sits on each ladder (fraction of the way up, in tens)")
foot={"as weighed (PDG)":([2.16,1273,172570],[4.70,93.5,4183]),
      "all at the Z (Xing, Zhang, Zhou 2012 type values)":([1.27,619,171700],[2.90,55,2890]),
      "all at 1e12 GeV (same source)":([0.58,281,82600],[1.24,24.0,1160])}  # MeV
nu=nus(0.36); N=[le[i]**2/nu[i] for i in range(3)]
print(f"   neutrinos {f(nu):.4f}   leptons {f(le):.4f}   hidden partners {f(N):.4f}")
for nm,(u,d) in foot.items(): print(f"   {nm:52s} up {f(u):.4f}  down {f(d):.4f}  sum {f(u)+f(d):.4f}  average {(f(u)+f(d))/2:.4f}")
# measurement error on the sum (PDG as weighed): mu 2.16+-.07(+.49-.26 really), md 4.70+-.07, ms 93.5+-.8, mc 1273+-4.6, mb 4183+-7, mt 172570+-290
S=[]
for _ in range(200000):
    u=[rng.normal(2.16,0.07*3),rng.normal(1273,4.6),rng.normal(172570,290)]; d=[rng.normal(4.70,0.07),rng.normal(93.5,0.8),rng.normal(4183,7)]
    if u[0]>0: S.append(f(u)+f(d))
S=np.array(S); print(f"   sum with measurement errors (up quark error tripled, as the PDG range is lopsided): {S.mean():.4f} +- {S.std():.4f}  -> distance from 1: {(S.mean()-1)/S.std():+.1f} sigma")
# chance: five ladders, 10 pairs, test 'two middles add to 1 within 0.006'
fs=[f(nu),f(le),f(N),f(foot['as weighed (PDG)'][0][0:3]),f(foot['as weighed (PDG)'][1])]
cnt=0; T=200000
for _ in range(T):
    r=rng.uniform(0.3,0.7,5); ok=False
    for i in range(5):
        for j in range(i+1,5):
            if abs(r[i]+r[j]-1)<0.006: ok=True
    cnt+=ok
print(f"   chance that some pair of five ladders adds to 1 this closely: {100*cnt/T:.0f}%   (for the up/down pair alone: {100*0.012/0.4:.0f}%)")
print("\n(2) Which surface does the neutrino spiral pick? (four steps below neutrino 1, own rule, nothing adjusted)")
H0eV=lambda h: h*2.1332e-33/100*1.0   # hbar*H0 in eV for H0 in km/s/Mpc: 100 km/s/Mpc = 2.1332e-33 eV
def surface(m1,dm21=74.2,dm31=2510.):
    m=[m1,(m1**2+dm21)**.5,(m1**2+dm31)**.5]; L=lg(np.array(m)*1e-3); d1,d2=L[1]-L[0],L[2]-L[1]; k=d2/d1; Linf=L[2]+d2*k/(1-k); C=Linf-L[0]; return 10**(Linf-C*k**(-4))
for m1 in (0.30,0.34,0.355,0.36,0.365,0.38,0.42):
    s=surface(m1); print(f"   lightest neutrino {m1:.3f} meV -> surface {s:.3g} eV = H0 {s/2.1332e-35:.3g} km/s/Mpc")
from scipy.optimize import brentq
for h in (67.4,73.0): print(f"   H0 = {h}: needs lightest neutrino {brentq(lambda m: lg(surface(m))-lg(h*2.1332e-35),0.2,0.6):.4f} meV")
# error from the splittings (NuFIT: dm21 74.2+-2.1, dm31 2510+-27 in 1e-6 eV^2) at fixed m1 = 0.36
X=[lg(surface(0.36,rng.normal(74.2,2.1),rng.normal(2510,27))/2.1332e-35) for _ in range(20000)]
print(f"   at 0.36 meV exactly, splitting errors alone move H0 by a factor of 10^{np.std(X):.2f} = x{10**np.std(X):.1f}")
print("   -> the spiral cannot tell 67 from 73; it fixes the surface only to within a factor of a few.")
print("\n(3) Both Hubble values at once")
print(f"   ratio 73.0/67.4 = {73.0/67.4:.4f} (8.3% apart); the early value assumes the sound-horizon ruler 147.1 Mpc; if both are right the ruler is {147.1*67.4/73.0:.1f} Mpc ({100*(67.4/73.0-1):.1f}%)")
print(f"   a shorter ruler needs faster early expansion: extra energy before atoms formed of about {100*((73.0/67.4)**2-1):.0f}% for a short time (rough: H^2 ~ energy)")
