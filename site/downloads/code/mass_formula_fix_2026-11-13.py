"""TRPG mass-formula fix (13 Nov 2026).
1) Count the adjustable numbers of the earlier mass formula. 2) Koide (1981) and Brannen (2006) as a zero/one-knob replacement in sqrt-mass space,
which is the S_ij = 1/2 ln(mi/mj) space (sqrt(mi/mj) = e^S). 3) Random-mass test of how special the result is.
Masses: PDG 2024 (leptons pole; quarks MS-bar u,d,s @2GeV, c,b @own mass, t pole)."""
import numpy as np, math
me,mmu,mtau=0.51099895,105.6583755,1776.93
def Q(a,b,c): return (a+b+c)/(math.sqrt(a)+math.sqrt(b)+math.sqrt(c))**2
print("Knobs of the earlier formula: E0 (fit to proton), P0 (soliton), '/10', anchors 17 & 40, v=174.1 GeV, O(1) coeffs 0.218,0.415,0.439,1.71,2.66,0.9925 -> >=11 choices for 9 masses")
print("Koide Q leptons = %.7f (2/3 = %.7f)"%(Q(me,mmu,mtau),2/3))
# predict tau from e, mu with Q=2/3: solve for x=sqrt(mtau)
a,b=math.sqrt(me),math.sqrt(mmu); s=a+b
# (a^2+b^2+x^2) = 2/3 (s+x)^2 -> x^2/3 - 4/3 s x + (a^2+b^2) - 2/3 s^2 = 0
A_,B_,C_=1/3,-4/3*s,(a*a+b*b)-2/3*s*s
x=(-B_+math.sqrt(B_*B_-4*A_*C_))/(2*A_); print("tau predicted from e,mu: %.3f MeV vs %.2f measured (%.4f%%)"%(x*x,mtau,(x*x/mtau-1)*100))
# Brannen form sqrt(m_k)=M(1+sqrt2 cos(d+2pi k/3)); fit d from data, compare to 2/9
r=np.sqrt([me,mmu,mtau]); M=r.mean()
d=[math.acos((ri/M-1)/math.sqrt(2)) for ri in r]
print("Brannen scale M^2 = %.4f MeV; angle data:"%(M*M), [round(v,5) for v in d])
for dd in [2/9]:
    pred=[ (M*(1+math.sqrt(2)*math.cos(dd+2*math.pi*k/3)))**2 for k in range(3)]
    print("delta=2/9, one scale -> masses",[round(p,4) for p in sorted(pred)])
# Random test: log-uniform triples over 1e-1..1e4 MeV, chance |Q-2/3| <= |Q_lep-2/3|
rng=np.random.default_rng(1); N=2_000_000
m=10**rng.uniform(-1,4,(N,3)); q=m.sum(1)/np.sqrt(m).sum(1)**2
tol=abs(Q(me,mmu,mtau)-2/3); print("random chance of matching 2/3 this closely: %.2e"%np.mean(abs(q-2/3)<=tol))
quarks={"u":2.16,"d":4.70,"s":93.5,"c":1273,"b":4183,"t":172570}
for trip in [("u","c","t"),("d","s","b"),("s","c","b"),("c","b","t")]:
    print(trip,"Q=%.4f"%Q(*[quarks[k] for k in trip]))
# One-input version: fix delta=2/9, take scale only from the tau, predict e and mu
d=2/9; Mt=math.sqrt(mtau)/(1+math.sqrt(2)*math.cos(d))
pm=(Mt*(1+math.sqrt(2)*math.cos(d+4*math.pi/3)))**2; pe=(Mt*(1+math.sqrt(2)*math.cos(d+2*math.pi/3)))**2
print("from tau alone: e=%.5f (%.3f%%), mu=%.4f (%.4f%%)"%(pe,(pe/me-1)*100,pm,(pm/mmu-1)*100))
print("Brannen scale x3 = %.2f MeV vs proton 938.27 (%.2f%%)"%(3*Mt*Mt,(3*Mt*Mt/938.272-1)*100))
