"""One neutrino, three perspectives; and does neutrino mass depend on temperature/activity? (14 Nov 2026)
Data: NuFIT 6.0 (2024, normal ordering); record M2 masses; Borexino; Super-K; KamLAND; Daya Bay; T2K/NOvA; IceCube DeepCore.
Credit: Pontecorvo (1957, 1967); Maki, Nakagawa & Sakata (1962); Wolfenstein (1978); Mikheyev & Smirnov (1985);
Harrison, Perkins & Scott (2002, tribimaximal); Lam, Albright & Rodejohann (TM1); Koide; Brannen; Notzold & Raffelt (1988)."""
import math, numpy as np
s12, s13, s23 = 0.307, 0.02195, 0.470   # s23 also checked at 0.561
dm21, dm31 = 7.49e-5, 2.534e-3         # eV^2
m_rec = [0.00036, 0.0087, 0.0501]

print("="*78); print("A. Are the three neutrino masses one object seen three ways (no real splits)?"); print("="*78)
print(f"Record masses {m_rec}: dm21 = {m_rec[1]**2-m_rec[0]**2:.3e} (meas {dm21:.2e}), dm31 = {m_rec[2]**2-m_rec[0]**2:.3e} (meas {dm31:.3e}), sum {sum(m_rec):.4f} eV")
for E, L, name in [(3e-3, 53, "JUNO/reactor 3 MeV, 53 km"), (3e-3, 180, "KamLAND 3 MeV, 180 km"), (0.6, 295, "T2K 0.6 GeV, 295 km"), (25.0, 12700, "IceCube up-going 25 GeV")]:
    ph31 = 1.267*dm31*L/E; ph21 = 1.267*dm21*L/E     # L km, E GeV
    print(f"  {name:28s} oscillation phase (31) {ph31:7.2f} rad, (21) {ph21:6.3f} rad")
print("If the three masses were equal (one mass, three viewpoints only), every phase above is 0: no oscillation.")
print("Oscillation is seen, with the L/E dip (Super-K 2004, KamLAND 2008). So the three mass values are real and different.")
print("What CAN be one object: the neutrino is one field with a 3-state inside; flavours are views of it (this is standard).")

print("\n"+"="*78); print("B. Three perspectives a third of a turn apart"); print("="*78)
w = np.exp(2j*math.pi/3)
F = np.array([[w**(j*k) for k in range(3)] for j in range(3)])/math.sqrt(3)
print("Pure one-third-turn views: every |U|^2 = 1/3 -> sin^2 th13 = 0.333 vs 0.022 measured: EXCLUDED.")
TBM = np.array([[math.sqrt(2/3), math.sqrt(1/3), 0], [-math.sqrt(1/6), math.sqrt(1/3), -math.sqrt(1/2)], [-math.sqrt(1/6), math.sqrt(1/3), math.sqrt(1/2)]])
print("Tribimaximal (HPS 2002) = third-turn views with one pair merged: sin^2 th13 = 0: EXCLUDED (Daya Bay/RENO 2012).")
def angles(U):
    a13 = abs(U[0, 2])**2; a12 = abs(U[0, 1])**2/(1-a13); a23 = abs(U[1, 2])**2/(1-a13)
    J = np.imag(U[0, 0]*U[1, 1]*np.conj(U[0, 1])*np.conj(U[1, 0]))
    c12, c13, c23 = math.sqrt(1-a12), math.sqrt(1-a13), math.sqrt(1-a23)
    sd = J/(math.sqrt(a12)*c12*math.sqrt(a23)*c23*math.sqrt(a13)*c13**2)
    return a12, a13, a23, J, sd
# TM1: keep first column (2/3,1/6,1/6); misalign the other two views by angle th with phase ph
print("TM1: keep one view of the third-turn set, misalign the other two by one angle theta (and a phase):")
best = []
for ph in np.radians(np.arange(0, 360, 1)):
    # choose theta to reproduce measured s13
    lo, hi = 0.0, 0.5
    for _ in range(60):
        th = 0.5*(lo+hi)
        R = np.array([[1, 0, 0], [0, math.cos(th), math.sin(th)*np.exp(-1j*ph)], [0, -math.sin(th)*np.exp(1j*ph), math.cos(th)]])
        if abs((TBM@R)[0, 2])**2 < s13: lo = th
        else: hi = th
    a12, a13, a23, J, sd = angles(TBM@R)
    best.append((math.degrees(ph), math.degrees(th), a12, a23, J))
best = np.array(best)
print(f"  theta needed for th13: {best[0,1]:.2f} deg (one misalignment angle)")
print(f"  predicted sin^2 th12 = {best[0,2]:.4f} for every phase (measured 0.307 +/- 0.012): {(best[0,2]-0.307)/0.012:+.1f} sigma")
for target in (0.470, 0.561):
    i = np.argmin(abs(best[:, 3]-target)); J = best[i, 4]
    print(f"  if sin^2 th23 = {target}: needs phase {best[i,0]:.0f} deg; predicted Jarlskog J = {J:+.4f} (max possible about 0.034)")
print("  NuFIT 6.0 prefers delta_CP near 180-215 deg (J small, weakly measured). TM1 with s23 near 0.47-0.56 predicts |J| near maximal.")
print("  -> a real test for DUNE and Hyper-Kamiokande (delta_CP to about 10 deg).")

print("\n"+"="*78); print("C. Do the masses follow the lepton 'self = cross' rule (Koide/Brannen)?"); print("="*78)
def Q(ms, signs=(1, 1, 1)):
    r = [s*math.sqrt(m) for s, m in zip(signs, ms)]
    return sum(ms)/sum(r)**2
print(f"  record masses: Koide Q = {Q(m_rec):.4f} (charged leptons 0.6667; Q=2/3 means self = cross)")
print(f"  with the lightest root negative (Brannen 2006 style): Q = {Q(m_rec,(-1,1,1)):.4f}")
def solve(signs):
    out = []
    for m1 in np.logspace(-7, -1, 200001):
        m2 = math.sqrt(m1**2+dm21); m3 = math.sqrt(m1**2+dm31); out.append((m1, Q((m1, m2, m3), signs)-2/3))
    sol = []
    for (a, fa), (b, fb) in zip(out, out[1:]):
        if fa*fb < 0: sol.append(a)
    return sol
for sg in [(1, 1, 1), (-1, 1, 1)]:
    sols = solve(sg)
    for m1 in sols:
        m2 = math.sqrt(m1**2+dm21); m3 = math.sqrt(m1**2+dm31)
        print(f"  signs {sg}: Q=2/3 needs m1 = {m1:.5f}, m2 = {m2:.5f}, m3 = {m3:.5f}, sum = {m1+m2+m3:.4f} eV")
    if not sols: print(f"  signs {sg}: no solution with the measured splits")
print("  DESI 2025 + CMB: sum < ~0.064 eV.")

print("\n"+"="*78); print("D. Does neutrino mass depend on temperature or activity (energy)?"); print("="*78)
data = [("KamLAND reactor (~3 MeV, core ~600 K)", "dm21", 7.53e-5, 0.18e-5, 0.003),
        ("Solar (0.1-15 MeV, core 15 million K)", "dm21", 6.10e-5, 1.2e-5, 0.008),
        ("Daya Bay reactor (~4 MeV)", "dm31", 2.519e-3, 0.06e-3, 0.004),   # dm32 ~2.47 +/- converted to dm31
        ("T2K/NOvA accelerator (0.6-2 GeV)", "dm31", 2.52e-3, 0.04e-3, 1.5),
        ("Super-K atmospheric (0.1-100 GeV)", "dm31", 2.50e-3, 0.10e-3, 5.0),
        ("IceCube DeepCore (6-60 GeV)", "dm31", 2.48e-3, 0.06e-3, 25.0)]
for nm, which, v, e, E in data: print(f"  {nm:40s} {which} = {v:.3e} +/- {e:.1e} eV^2  (typical E {E} GeV)")
d31 = [(v, e, E) for _, w_, v, e, E in data if w_ == "dm31"]
vals = np.array([d[0] for d in d31]); errs = np.array([d[1] for d in d31]); Es = np.array([d[2] for d in d31])
x = np.log(Es); y = np.log(vals); wt = (vals/errs)**2
A = np.vstack([np.ones_like(x), x]).T; W = np.diag(wt)
cov = np.linalg.inv(A.T@W@A); p = cov@A.T@W@y
print(f"  fit dm31 proportional to E^p across {Es.min()*1000:.0f} MeV - {Es.max():.0f} GeV: p = {p[1]:+.4f} +/- {math.sqrt(cov[1,1]):.4f}")
print(f"  -> mass^2 changes by at most {100*(math.exp(2*math.sqrt(cov[1,1])*math.log(Es.max()/Es.min()))-1):.1f}% over a factor {Es.max()/Es.min():.0f} in energy (2 sigma)")
print("  The solar value sits about 1.2 sigma below KamLAND: consistent (it was ~2 sigma before 2020).")
print("Verdict: the vacuum mass splits do NOT depend on source temperature or neutrino energy.")
print("\nWhat DOES change with surroundings: the MSW matter effect (Wolfenstein 1978; Mikheyev & Smirnov 1985).")
GF = 1.1663787e-5; hbarc = 1.97327e-14  # GeV^-2, GeV*cm
for ne_cm3, place in [(6.0e25, "solar core (rho ~150 g/cm3)"), (1.5e24, "Earth's mantle"), (5e22*0+3e23, "outer Sun")]:
    V = math.sqrt(2)*GF*ne_cm3*hbarc**3*1e9   # eV
    E_res = dm21*math.cos(2*math.asin(math.sqrt(s12)))/(2*V)/1e6   # MeV
    print(f"  {place:28s}: extra potential V = {V:.2e} eV; solar-split resonance energy = {E_res:.1f} MeV")
c13s = (1-s13)**2; s2_12 = 4*s12*(1-s12)
print(f"  predicted survival low energy (vacuum-like): {c13s*(1-0.5*s2_12):.2f};  high energy (matter-dominated): {c13s*s12+s13**2:.2f}")
print("  Borexino measured: pp (low) 0.57 +/- 0.09; 8B (high) about 0.37 +/- 0.08 -> the environment effect is seen.")
print("  Super-K day/night (Earth matter at night): -3.3% +/- 1.1% (2016).")
print("  Hot plasma (early universe, supernova): neutrinos gain a thermal effective mass ~ G_F T^4 terms (Notzold & Raffelt 1988).")
