# Do whole waves fit between the surface and the grain, and does a value come out?
import numpy as np
lg=np.log10
le=[0.51099895e6,105.6584e6,1776.93e6]; up=[2.16e6,1273e6,172570e6]; dn=[4.70e6,93.5e6,4183e6]
def nus(m1): return [m1*1e-3,((m1*1e-3)**2+7.42e-5)**.5,((m1*1e-3)**2+2.51e-3)**.5]
GRAIN,MPL,H0,MMG,DD,DE=8.63e27,1.2209e28,1.44e-33,2.41e17,4.62e22,2.24e-3
LS,LG=lg(H0),lg(GRAIN); gm=lambda m: 10**np.mean(lg(np.array(m)))
nu=nus(0.36); N=[le[i]**2/nu[i] for i in range(3)]
print("(1) Stretched and squished, literally: wavelength = hc/E")
for n,E in (("surface",H0),("neutrino 1",nu[0]),("electron",le[0]),("top",up[2]),("hidden 3",N[2]),("grain",GRAIN)): print(f"   {n:10s} {E:.3g} eV  wavelength {1.2398e-6/E:.3g} m")
print("\n(2) Whole turns between the ends: surface to grain = %.2f decades"%(LG-LS))
for n in (3,4,5,6,7,8):
    J=(LG-LS)/n; print(f"   {n} turns: one turn = x{10**J:.3g}; where the ladder middles fall, in turns from the surface:", {k:round((lg(gm(m))-LS)/J,3) for k,m in (("nu",nu),("lep",le),("down",dn),("up",up),("hidden",N))})
J=(LG-LS)/6
print("\n(3) Six turns: the levels that come out, against what sits there")
lev=[10**(LS+k*J) for k in range(7)]
for k,(what,val) in zip(range(7),[("the surface",H0),("nothing known",None),("nothing known",None),("neutrino middle (m1 0.36) / dark energy",(gm(nu),DE)),("lepton middle",gm(le)),("hidden middle (m1 0.36) / MMG",(gm(N),MMG)),("the grain",GRAIN)]):
    print(f"   level {k}: {lev[k]:.3g} eV   {what}: {val}")
print("   lepton middle / level 4 = %.3f ; in turns from the surface: %.4f"%(gm(le)/lev[4],(lg(gm(le))-LS)/J))
print("   = grain^(2/3) x surface^(1/3)  (Weinberg 1972 form (H0 Mpl^2)^(1/3): with full Planck energy %.3g eV)"%((H0*MPL**2)**(1/3)))
print("   sensitivity: H0 67-73 km/s/Mpc ->",[round((lg(gm(le))-lg(h))/((LG-lg(h))/6),4) for h in (1.43e-33,1.56e-33)],"; grain x2 / x0.5 ->",[round((lg(gm(le))-LS)/((lg(g)-LS)/6),4) for g in (GRAIN*2,GRAIN/2)])
print("   neutrino middle on level 3 needs m1 = ",end="")
from scipy.optimize import brentq
print("%.3f meV (the spiral landings need 0.36)"%brentq(lambda m: lg(gm(nus(m)))-lg(lev[3]),0.01,1))
print("   the lepton triple predicted from level 4 with the 45 deg rule and the twist: scale only; e,mu,tau all x%.3f"%(lev[4]/gm(le)))
# chance: a middle within this distance of some level, for 3..8 turns tried
x=(lg(gm(le))-LS)/J; d=abs(x-round(x)); print("   lepton middle is %.4f turn from a whole turn; chance for one fixed count %.1f%%; allowing 3 to 8 turns to be tried ~%.0f%%"%(d,2*d*100,min(100,2*d*100*6)))
print("\n(4) Each ladder's own spiral: does it meet an end in a whole number of steps?")
for n,m in (("neutrinos",nu),("leptons",le),("hidden",N),("up",up),("down",dn)):
    L=lg(np.array(m)); d1,d2=L[1]-L[0],L[2]-L[1]; k=d2/d1; Li=L[2]+d2*k/(1-k); C=Li-L[0]
    end=LS if k<1 else LG; t=np.log((Li-end)/C)/np.log(k); print(f"   {n:10s} meets the {'surface' if k<1 else 'grain  '} at step {t:+.2f}")
print("\n(5) If every wave had to squish AT THE GRAIN (distance to the grain shrinking by a fixed factor): middle rung predicted from the two ends")
for n,m in (("neutrinos",nu),("leptons",le),("hidden",N),("up",up),("down",dn)):
    L=lg(np.array(m)); p=LG-np.sqrt((LG-L[0])*(LG-L[2])); print(f"   {n:10s} predicted middle {10**p:.3g} eV, actual {m[1]:.3g} eV (x{m[1]/10**p:.2f})")
