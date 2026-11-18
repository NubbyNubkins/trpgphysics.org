"""TRPG master action, step M2: the mass couplings Y as a relational 3-generation matrix (13 Nov 2026).
sqrt-mass matrix  S = M0 * [[1, r w, r w*],[r w*, 1, r w],[r w, r w*, 1]],  w = e^{i delta}  (cyclic / circulant)
eigenvalues M0(1 + 2r cos(delta + 2 pi k/3)); Koide Q=2/3  <=>  r = 1/sqrt2.
Refs: Koide 1981/1983; Foot 1994; Brannen 2006. Neutrino data: NuFIT 6.0 (2024). Cosmology: DESI DR2 (2025)."""
import numpy as np, math
me,mmu,mtau=0.51099895,105.6583755,1776.93
r=1/math.sqrt(2); d=2/9
M0=(math.sqrt(me)+math.sqrt(mmu)+math.sqrt(mtau))/3
w=np.exp(1j*d); S=M0*np.array([[1,r*w,r*np.conj(w)],[r*np.conj(w),1,r*w],[r*w,r*np.conj(w),1]])
ev=np.sort(np.linalg.eigvalsh(S)); print("matrix -> masses",np.round(ev**2,4),"measured",[me,mmu,mtau])
# 'equal self vs relational' check: |common part|^2 vs |relational part|^2 of the sqrt-mass vector
v=np.sqrt([me,mmu,mtau]); c=v.mean()*np.ones(3); p=v-c
print("self (common) part^2 = %.3f, relational part^2 = %.3f  (Koide = equal)"%(c@c,p@p))
# Neutrino test: impose Q=2/3 with measured splittings (normal ordering), allow one negative sqrt (as for Brannen)
dm21,dm31=7.49e-5,2.513e-3   # eV^2, NuFIT 6.0 normal ordering
from scipy.optimize import brentq
def Q(m,s): rt=np.array(s)*np.sqrt(m); return m.sum()/rt.sum()**2
for s in [(1,1,1),(-1,1,1)]:
    f=lambda m1: Q(np.array([m1,math.sqrt(m1**2+dm21),math.sqrt(m1**2+dm31)]),s)-2/3
    grid=np.logspace(-6,0,4000); vals=[f(x) for x in grid]
    roots=[brentq(f,grid[i],grid[i+1]) for i in range(len(grid)-1) if np.sign(vals[i])!=np.sign(vals[i+1])]
    for m1 in roots:
        ms=[m1,math.sqrt(m1**2+dm21),math.sqrt(m1**2+dm31)]
        print("signs",s,"-> m1,m2,m3 = %.5f %.5f %.5f eV, sum = %.4f eV"%(*ms,sum(ms)))
    if not roots: print("signs",s,"-> no solution")
# Inverted ordering check: m3 lightest, m1 = sqrt(m3^2+|dm32|+...), use dm32=-2.484e-3
dm32=-2.484e-3
for s in [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)]:
    def f(m3,s=s):
        m2=math.sqrt(m3**2-dm32); m1=math.sqrt(m2**2-dm21); return Q(np.array([m1,m2,m3]),s)-2/3
    grid=np.logspace(-6,0,4000); vals=[f(x) for x in grid]
    roots=[brentq(f,grid[i],grid[i+1]) for i in range(len(grid)-1) if np.sign(vals[i])!=np.sign(vals[i+1])]
    print("inverted signs",s,"roots (m3 eV):",[round(x,5) for x in roots])
