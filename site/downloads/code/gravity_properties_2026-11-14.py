"""Properties 1 (long range) and 2 (acts on all energy equally) for TRPG's resistance field. 14 Nov 2026."""
hc_eVm=1.23984198e-6
for name,m_eV in [("containment field ripple P (M1, glueball)",1.72e9),("graviton upper limit (LIGO-Virgo-KAGRA GWTC-3, 2021)",1.27e-23)]:
    lam=hc_eVm/m_eV; print("%-55s mass %.2e eV -> reach ~ %.2e m"%(name,m_eV,lam))
print("reach of P = %.2f fm  -> P cannot be gravity (too short by > 10^30)"%(hc_eVm/1.72e9*1e15))
print("GW speed vs light (GW170817): equal to within ~3e-15 -> consistent with massless")
print("MICROSCOPE (2022): different materials fall alike to 1.5e-15 -> every TRPG layer must feel the same resistance to that level")
