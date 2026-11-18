# Temperature (not time alone) as the driver; the majoron as early resistance; the start level; ages against the Hubble values
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
ln=np.log; pi=np.pi
MPL=1.2209e19  # GeV
def H(T,g): return 1.66*np.sqrt(g)*T*T/MPL          # GeV, radiation era, T in GeV
def age(T,g): return 2.42/np.sqrt(g)/(T*1e3)**2      # s
eta=6.1e-10; zeta3=1.2020569
print("(A) What stops each creation: temperature, plus two other things")
# 1 neutrinos let go: weak rate G_F^2 T^5 against the growth rate of space
GF=1.1664e-5; Tnu=(1.66*np.sqrt(10.75)/(GF**2*MPL))**(1/3)
print(f"  neutrinos let go when their meeting rate falls below the growth rate of space: {Tnu*1e3:.2f} MeV (their own energy is 0.00000036 to 0.00005 MeV: not set by their weight at all)")
# 2 first nuclei: deuterium survives when light bits energetic enough to break it are fewer than matter bits
B=2.2246e-3; mp=0.93827
fD=lambda T: ln(eta)+1.5*ln(T/mp)+B/T+ln(12*zeta3/np.sqrt(pi))  # ~ n_D/n_p of order 1 (Saha, order-one factors kept roughly)
TD=brentq(fD,0.03e-3,0.5e-3); print(f"  first nuclei: binding 2.22 MeV, but they hold only below {TD*1e3:.3f} MeV = binding / {B/TD:.0f}")
# 3 atoms: Saha, half the electrons bound
me=0.51099895e6; BH=13.6057
def saha(T,x=0.5,et=eta):   # T in eV ; x^2/(1-x) = (1/n_b) (me T/2pi)^1.5 exp(-B/T), n_b = eta * 2 zeta3/pi^2 T^3 (hydrogen fraction 0.76 of baryons by number ~0.92; use 0.88 e per baryon avg -> keep 1)
    nb=et*2*zeta3/pi**2*T**3; return ln((me*T/(2*pi))**1.5/nb)-BH/T-ln(x*x/(1-x))
Trec=brentq(saha,0.1,1.0); print(f"  atoms: binding 13.6 eV, but half the electrons are held only below {Trec:.3f} eV = binding / {BH/Trec:.0f} (light finally set free a little later, 0.26 eV)")
print(f"  both delays have ONE cause: there are {1/eta*1:.2e} bits of light for every bit of matter (1 / the matter excess).")
for e in (6.1e-11,6.1e-10,6.1e-9):
    td=brentq(lambda T: ln(e)+1.5*ln(T/mp)+B/T+ln(12*zeta3/np.sqrt(pi)),0.03e-3,0.5e-3); tr=brentq(lambda T: saha(T,0.5,e),0.1,1.0)
    print(f"     if the matter excess were {e:.1e}: nuclei at {td*1e3:.3f} MeV ({age(td,3.36)/60:.1f} min), atoms at {tr:.3f} eV")
print("  quarks lock into protons at 155 MeV and the surface's resistance switches on at 160 GeV: these two ARE set by temperature alone (the balance point melts above it).")
print("\n(B) The majoron as the early resistance for the Hubble values")
f=6.3e10; mnu3=0.0501e-9  # GeV
lam=mnu3/f; print(f"  main line: settling value f = {f:.2g} GeV -> the majoron's grip on neutrino 3 = weight / f = {lam:.2e}")
for m in (0.1,1.0):   # eV
    G=lam**2*(m*1e-9)/(16*pi); Hh=H(m*1e-9/3,3.36); print(f"  majoron of {m} eV: meeting rate / growth rate of space at its own temperature = {G/Hh:.1e}   (needs about 1)")
print(f"  needed grip about 1e-13 (Escudero & Witte 2020) -> f would have to be {mnu3/1e-13:.0f} GeV, not 6.3e10: off by a factor {f/(mnu3/1e-13):.1e}")
print(f"  and its weight: 8.8 keV / angle^4 with the largest angle pi gives at least {8800/pi**4:.0f} eV, not 0.1 to 1 eV.  -> FAIL twice.")
print("\n(C) Ages against the two Hubble values (flat, matter + constant dark energy, same early universe: matter x h^2 = 0.1430)")
def t0(h):
    Om=0.1430/h**2; Or=4.18e-5/h**2; OL=1-Om-Or; return quad(lambda a: 1/(a*np.sqrt(Or/a**4+Om/a**3+OL)),0,1)[0]*977.8/(100*h)
for h in (0.674,0.70,0.73): print(f"  H0 = {100*h:.1f}: age today {t0(h):.2f} billion years")
print("  oldest star clusters: 13.5 +- 0.3 (stat) +- 0.5 (sys) billion years after adding the ~0.2 before stars (Valcin et al. 2021)")
print("\n(D) The start level: where on the cone the light began")
N=[7.25e5,1.29e9,6.30e10]; MMG=2.41e8
for th in (0.5,1.0,1.41,2.0,3.0):
    m=8.8e-6/th**4; Tosc=np.sqrt(m*MPL/(3*1.66*np.sqrt(106.75))); print(f"  starting angle {th}: majoron {m*1e6:.3g} keV begins to swing when the light is at {Tosc:.3g} GeV (age {age(Tosc,106.75):.1e} s)")
print(f"  window: above the swing level, below the MMG level {MMG:.3g} GeV (age {age(MMG,106.75):.1e} s) so the MMG are never made")
# hidden partner 1 is inside the window: thermal. Does it erase the excess made earlier?
K=0.36/1.07; print(f"  hidden partner 1 ({N[0]:.3g} GeV) is inside the window, so it IS made by the light. Erasing strength K = 0.36 meV / 1.07 meV = {K:.2f}: it leaves {100*np.exp(-3*pi*K/8):.0f}% of an earlier excess")
v=174.1; eps=3*N[0]*0.0501e-9/(16*pi*v*v); print(f"  could hidden partner 1 make the excess itself? largest lopsidedness {eps:.1e} (Davidson & Ibarra 2002) -> excess at most {eps*0.96e-2:.1e} against 6.1e-10 measured: short by x{6.1e-10/(eps*0.96e-2):.0f}: NO")
for i,n in enumerate(N):
    l=n/f; G=l**4*n/(8*pi*16); print(f"  majorons heated by hidden partner {i+1}? rate/growth at its level = {G/H(n,106.75):.1e} {'(would heat them: must stay unmade)' if G/H(n,106.75)>1 else '(no)'}")
