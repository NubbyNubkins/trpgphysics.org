"""Neutron stars with the containment pressure B set by the proton: finer scan of the quark pairing gap.

Uses the equation of state and solver of neutron_stars_no_knobs_2026-11-14.py (same folder): perturbative-QCD
quark matter with a4 = 1 - 2 alpha_s/pi (Freedman and McLerran 1977; Fraga, Pisarski and Schaffner-Bielich 2001),
a colour-superconducting pairing gap (Alford, Rajagopal and Wilczek 1998; Rapp, Schafer, Shuryak and Velkovsky 1998;
Alford, Braby, Paris and Reddy 2005), the TOV equations and Hinderer's (2008) tidal deformability.
Limits: maximum mass at least 2.0 solar masses (PSR J0740+6620); Lambda(1.4) below about 580 (GW170817).
"""
import importlib.util, os
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("nk", os.path.join(here, "neutron_stars_no_knobs_2026-11-14.py"))
nk = importlib.util.module_from_spec(spec); spec.loader.exec_module(nk)

if __name__ == "__main__":
    print(" B (MeV/fm3)  gap (MeV)   Mmax (Msun)   R1.4 (km)   Lambda1.4   2.0 test   tidal test")
    for B, gaps in ((57.1, (70, 75, 80, 85, 90)), (60.1, (70, 80, 85, 90, 95, 100)), (63.7, (80, 90, 95, 105, 110))):
        for gap in gaps:
            M, L, R, a4, es = nk.sequence(B, 3.0, float(gap))
            print(f"   {B:5.1f}       {gap:4d}        {M:.2f}         {R:5.2f}       {L:5.0f}      "
                  f"{'PASS' if M >= 2.0 else 'fail'}       {'PASS' if L < 580 else 'fail'}", flush=True)
