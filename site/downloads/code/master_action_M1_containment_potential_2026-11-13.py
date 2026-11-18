"""TRPG master action, step M1: pin the containment ('persistence'/space-resistance) potential V(P) with no tuning.
V(P) = (lam/4) v^4 (P^2-1)^2 ; P=1 empty space (full resistance), P=0 inside a proton (resistance removed).
Inputs: B from neutron-star + tidal window (K1,K2) and the lightest glueball mass from lattice QCD
(Morningstar & Peardon 1999; Chen et al. 2006: 0++ ~ 1.71-1.73 GeV). Model form: Friedberg & Lee (1977)."""
import math
hbarc=197.3269804; MeV4_per_MeVfm3=hbarc**3
for B in (57.1,60.1,63.7,62,70):   # proton-set band (proton_sets_B) and the old neutron-star window
    for mg in (1710.,1730.):
        Bn=B*MeV4_per_MeVfm3; lam=mg**4/(16*Bn); v=math.sqrt(mg**2/(2*lam))
        wall=hbarc/mg
        print(f"B={B} MeV/fm3  glueball={mg:.0f} MeV -> lam={lam:7.0f}  v={v:5.1f} MeV  wall thickness ~{wall:.3f} fm")
# Thin wall (0.11 fm) << core radius (~0.95 fm bag, 0.7 fm rms) -> MIT bag limit is valid, so K1-K3 carry over unchanged.
# Compare with the QCD vacuum energy from the gluon condensate (Shifman, Vainshtein & Zakharov 1979):
G2=0.012 # GeV^4, <alpha_s G^2/pi>
eps=(9/32)*G2*1e12/MeV4_per_MeVfm3
print(f"SVZ vacuum energy ~ {eps:.0f} MeV/fm3, about {eps/66:.1f}x the bag B (a known, unresolved gap)")
