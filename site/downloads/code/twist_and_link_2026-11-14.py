"""The 2/9 twist and 1/sqrt2 link in the lepton matrix (M2). 14 Nov 2026. PDG 2024 masses and errors."""
import math, numpy as np
def delta(me,mmu,mtau):
    r=np.sqrt([me,mmu,mtau]); M=r.mean(); return math.acos((r[2]/M-1)/math.sqrt(2))
d0=delta(0.51099895,105.6583755,1776.93)
rng=np.random.default_rng(0)
ds=[delta(0.51099895,105.6583755+rng.normal(0,2.3e-6),1776.93+rng.normal(0,0.09)) for _ in range(20000)]
print("twist from data: %.6f +/- %.6f rad ; 2/9 = %.6f ; difference = %.1f sigma"%(d0,np.std(ds),2/9,(d0-2/9)/np.std(ds)))
Q=2/3; print("2/9 = Q/3 = (Koide ratio)/(number of generations) = %.6f"%(Q/3))
print("link 1/sqrt2: Koide Q=2/3 <=> shared part^2 = relational part^2 (shown in M2); r=1/sqrt2 is that balance")
