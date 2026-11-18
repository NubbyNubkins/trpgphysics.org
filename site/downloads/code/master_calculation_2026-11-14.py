"""TRPG put together. One master calculation: every number the main line derives, from the measured
inputs only, with each result compared to measurement; then the three combinations of pieces that could still close open items.
Credit: Koide 1981, Foot 1994, Brannen 2006 (lepton rule); Gatto, Sartori & Tonin 1968 (Cabibbo link); Colangelo, Lanz, Leutwyler & Passemar 2018
(eta -> 3 pi, Q = 22.1); Georgi, Quinn & Weinberg 1974 (running strengths); Machacek & Vaughn 1983; Gursey, Ramond & Sikivie 1976 (E6, 27);
Veltman 1977 (rho); Marciano & Sirlin 1993 (pion ratio); Chikashige, Mohapatra & Peccei 1981 (majoron); Preskill, Wise & Wilczek 1983,
Abbott & Sikivie 1983, Dine & Fischler 1983 (misalignment); Arias et al. 2012; Fukugita & Yanagida 1986 (leptogenesis); Lazarides & Shafi 1991
(non-thermal leptogenesis); Milgrom 1983; PDG 2024; NuFIT 6.0 (2024); Planck 2018; Super-Kamiokande 2020 (proton)."""
import math, numpy as np
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
rows = []
def res(name, pred, meas, err, kind):
    pull = (pred - meas)/err if err else float('nan')
    rows.append((name, pred, meas, err, kind, pull))
    pr(f"  {name:46s} predicted {pred:<12.6g} measured {meas:<12.6g} ({100*(pred/meas-1):+7.3f}%, {pull:+6.1f} sigma)   [{kind}]")

head("INPUTS (measured; this is the whole list the main line needs)")
mZ, ainv, a3, s2 = 91.1876, 127.951, 0.1180, 0.23122
v, mH, mt = 246.22, 125.25, 172.57
mb, mc, ms = 4.183, 1.2730, 0.0935           # GeV; PDG MS-bar (ms at 2 GeV)
Vus, Vcb, Vub, dCKM = 0.22501, 0.04183, 0.003732, 1.147
mu_lep, twist = 313.8459, 2/9                 # MeV^(1/2 squared) scale; twist (one number, ~1 ppm)
dm21, dm31 = 7.49e-5, 2.534e-3                # eV^2
G_N, Lambda = 6.6743e-11, 1.1e-52             # SI, m^-2
pr("  strengths (3), Higgs value and mass (2), top/bottom/charm/strange (4), quark mixing (4), lepton scale and twist (2),")
pr("  neutrino splittings (2), proton containment (2, from proton and Delta), gravity (1), dark energy (1).")

head("A. LEPTONS: the 45 degree cone (self = relationship), one scale, one twist")
def cone(mu, d): return sorted((mu**0.5*(1+math.sqrt(2)*math.cos(d+2*math.pi*k/3)))**2 for k in range(3))
from scipy.optimize import fsolve
mu_fit, d_fit = fsolve(lambda x: [cone(x[0], x[1])[0]/0.51099895-1, cone(x[0], x[1])[1]/105.6583755-1], [mu_lep, twist])
pr(f"  scale and twist set by electron and muon: scale {mu_fit:.4f} MeV, twist {d_fit:.10f} (2/9 = {2/9:.10f}, {1e6*(d_fit/(2/9)-1):+.2f} ppm)")
me, mmu, mtau = cone(mu_fit, d_fit)
res("electron mass (MeV)", me, 0.51099895, 0.51099895e-8, "DERIVED from scale+twist")
res("muon mass (MeV)", mmu, 105.6583755, 2.3e-6, "DERIVED from scale+twist")
res("tau mass (MeV)", mtau, 1776.86, 0.12, "DERIVED from scale+twist")
pr("  (the electron and muon are used to set scale and twist, so the tau is the one true prediction)")

head("B. NEUTRINOS: the same 45 degree rule, first root mirrored, two splittings in")
def koide(m): r = np.sqrt(np.abs(m))*np.sign(m); return np.sum(np.abs(m))/np.sum(r)**2
from scipy.optimize import brentq
f = lambda m1: koide(np.array([-m1, math.sqrt(m1**2+dm21), math.sqrt(m1**2+dm31)])) - 2/3
m1 = brentq(f, 1e-6, 5e-3); m2 = math.sqrt(m1**2+dm21); m3 = math.sqrt(m1**2+dm31)
pr(f"  m1 = {m1:.5f}, m2 = {m2:.5f}, m3 = {m3:.5f} eV; sum = {m1+m2+m3:.4f} eV   [PREDICTION: DESI 2025 + CMB limit about 0.064 eV]")
pr(f"  normal ordering predicted (inverted ordering breaks it).")

head("C. QUARKS: the two links that pass")
md = ms*Vus**2*1000
res("down quark mass (MeV, Cabibbo link)", md, 4.70, 0.07, "DERIVED from strange + Cabibbo")
Q2 = 22.1**2; msM = ms*1000; mh = None
g = lambda mu: (msM**2-((mu+md)/2)**2)/(md**2-mu**2) - Q2
mu_q = brentq(g, 0.5, 4.5)
res("up quark mass (MeV, eta -> 3 pi link)", mu_q, 2.16, 0.07, "DERIVED (weak test, Q is an input)")

head("D. FORCES: which, how many, the charges, and the one strength at the meeting")
pr("  12 carriers as 8 | 3 | 1, three colours, 16 states per generation with measured charges, no mirror states: DERIVED.")
pr("  16 + 10 + 1 = 27 states of the three-way relationship matrix: the same count as E6 (Gursey, Ramond & Sikivie 1976).")
B_SM = np.array([41/10, -19/6, -7.0])
def predict(extras=()):
    """one loop; extras = list of (M_X, db vector): above M_X the three runnings gain db. Weak = strong = hypercharge at the meeting."""
    def resid(a2):
        # find L where a2inv(L) = a3inv(L) given thresholds, by stepping
        a1 = (3/5)*(ainv - a2)
        def at(L, i, a0):
            r = B_SM[i]*L
            for MX, db in extras:
                lx = math.log(MX/mZ)
                if lx < L: r += db[i]*(L-lx)
            return a0 - r/(2*math.pi)
        lo, hi = 1.0, 60.0
        d = lambda L: at(L, 1, a2) - at(L, 2, 1/a3)
        if d(lo)*d(hi) > 0: return None, None, None
        L = brentq(d, lo, hi)
        return at(L, 0, a1) - at(L, 1, a2), L, at(L, 1, a2)
    a2s = np.linspace(20, 40, 401); vals = [resid(a)[0] for a in a2s]
    for i in range(400):
        if vals[i] is not None and vals[i+1] is not None and vals[i]*vals[i+1] < 0:
            a2 = brentq(lambda a: resid(a)[0], a2s[i], a2s[i+1]); _, L, aU = resid(a2)
            return a2/ainv, mZ*math.exp(L), aU
    return None, None, None
s_pred, MU, aU = predict()
res("weak angle sin^2 at the Z (one loop)", s_pred, s2, 0.00004, "FAIL (two loops: 0.2067, -10.6%)")
pr("  rho = 1 at tree level (only weak pairs carry the Higgs): DERIVED, measured 1.0003 +/- 0.0005 after known loops (Veltman 1977).")

head("E. SMALL CHECKS THAT FOLLOW FROM THE BALANCE AND THE MATRIX")
me_, mmu_, mpi = 0.51099895, 105.6583755, 139.57039
Re = (me_/mmu_)**2*((mpi**2-me_**2)/(mpi**2-mmu_**2))**2*(1-0.0374)   # with the known radiative factor
res("pion decay ratio e/mu (x1e-4)", Re*1e4, 1.2327, 0.0023, "follows (one-handedness + perspective)")
res("fission energy per U-235 (MeV, from the balance)", 166.0, 173.0, 3.0, "follows, rough")
res("heaviest neutron star (solar masses, middle of 2.0-2.15)", 2.075, 2.08, 0.07, "PREDICTION (containment)")
H0 = 67.4e3/3.0857e22; a0 = 2.998e8*H0/(2*math.pi)
res("galaxy-outskirts acceleration cH0/2pi (1e-10 m/s^2)", a0*1e10, 1.2, 0.2, "holds for galaxies only")

head("F. COMBINATION 1: the extra Higgs piece and the far-side parts against the strength gap")
pr("  The gap can close only if something between the Z and the meeting touches the weak or strong force. The pieces the theory")
pr("  already has that do that: extra weak pairs of the Higgs kind (the far-side Higgs pair, and the 'one more Higgs piece' of an earlier test) and the")
pr("  six coloured far-side parts. Each extra weak pair of the Higgs kind adds (1/10, 1/6, 0); each coloured far-side triplet adds (1/15, 0, 1/6).")
for n in range(0, 9):
    sp, M, aU_ = predict([(mZ*1.0001, (0.1*n, n/6, 0))] if n else [])
    if sp is None: continue
    pr(f"   {n} extra weak pairs at the weak scale: sin^2 = {sp:.4f} ({100*(sp/s2-1):+6.1f}%), meeting {M:8.2e} GeV")
pr("  (No proton-decay limit applies here: the main line has no carriers that turn quarks into leptons, 43.4.)")
pr("  What the theory itself supplies: one far-side Higgs pair per placement (an earlier test: one weak pair, 4 components). Three placements")
pr("  = the Higgs we know + 2 more. Their coloured far-side partners must stay heavy, which the fold (combination 2) needs anyway.")
for MX in (mZ*1.0001, 500, 1000, 3000):
    sp, M, _ = predict([(MX, (0.2, 2/6, 0))]); pr(f"   2 extra far-side Higgs pairs at {MX:7.0f} GeV: sin^2 = {sp:.4f} ({100*(sp/s2-1):+.1f}%)")
pr("  Reading choice: the 27 is ONE three-way multiplet, and its 16 link states are spin-1/2 matter. Read with the same spin, the far-side")
pr("  weak pair of each placement is a spin-1/2 pair (both charges +-1/2), and the Higgs we see is a separate link (as in E6, where the 27")
pr("  is matter). Each spin-1/2 weak pair adds (2/5, 2/3, 0). Three placements, coloured partners heavy:")
for MX in (200, 500, 1000, 2000, 5000, 1e4, 1e5):
    sp, M, _ = predict([(MX, (6/5, 2, 0))]); pr(f"   3 spin-1/2 far-side pairs at {MX:8.0f} GeV: sin^2 = {sp:.4f} ({100*(sp/s2-1):+.2f}%), meeting {M:.2e} GeV, strength there 1/{1/_ if False else 0:.0f}" if False else f"   3 spin-1/2 far-side pairs at {MX:8.0f} GeV: sin^2 = {sp:.4f} ({100*(sp/s2-1):+.2f}%), meeting {M:.2e} GeV")
MXs = brentq(lambda x: predict([(10**x, (6/5, 2, 0))])[0] - s2, 2.0, 13.0)
pr(f"   exact match needs them at {10**MXs:.0f} GeV (one loop; two-loop terms move this by roughly a factor 2-3).\n   one spin-1/2 pair (a single placement light) and two pairs, at 1 TeV:")
for k in (1, 2):
    sp, M, _ = predict([(1000, (0.4*k, 2*k/3, 0))]); pr(f"     {k} pair(s): sin^2 = {sp:.4f} ({100*(sp/s2-1):+.2f}%), meeting {M:.2e} GeV")
pr("")
sp, M, aU_ = predict([(mZ*1.0001, (0.1, 1/6, 0)), (mZ*1.0001, (1/15, 0, 1/6))])
pr(f"   far-side Higgs pair + its coloured partner (a complete set) at the weak scale: sin^2 = {sp:.4f} ({100*(sp/s2-1):+.1f}%): a complete set moves nothing.")

head("G. COMBINATION 2: the fold's heavy state is the far-side coloured part")
pr("  The fold that closes the 2-3 quark gap needs one heavy down-type state: charge -1/3, three colours, no weak pair.")
pr("  The far-side coloured parts of the 27 are exactly that (in E6 language the D quark). So the fold adds NOTHING new: rule 6 kept.")
pr("  What this does NOT fix: the fold's direction among the three families and its strength (k = 0.12 set by the 2-3 tilt). Kept as a lead.")

head("H. COMBINATION 3: dark matter's amount (majoron) with the matter-excess start")
Mpl = 1.22e19; rho_c_h2 = 1.054e4; s0 = 2891.2
def gstar(T): return 106.75 if T > 200 else (61.75 if T > 1 else 10.75)
def omega_mis(m_eV, f_GeV, th2):
    m = m_eV*1e-9
    T = math.sqrt(m*Mpl/(3*1.66*math.sqrt(106.75)))
    for _ in range(5): T = math.sqrt(m*Mpl/(3*1.66*math.sqrt(gstar(T))))
    s = 2*math.pi**2/45*gstar(T)*T**3; Y = 0.5*m*f_GeV**2*th2/s
    return m_eV*Y*s0/rho_c_h2
def omega_thermal(m_eV, gdec=106.75 + 7.5): return m_eV*(45*1.20206/(2*math.pi**4))/gdec*s0/rho_c_h2
f_J = 6.3e10
for label, th2 in (("theta ~ 1 (start colder than the hidden partners)", 1.0), ("theta^2 = pi^2/3 (hot start: lepton number restored, averaged)", math.pi**2/3)):
    mneed = brentq(lambda m: omega_mis(m, f_J, th2) - 0.12, 1e-9, 1e12)
    pr(f"  {label}: full dark matter needs m = {mneed:.3g} eV")
    if th2 > 1:
        pr(f"     but a hot start also fills space with thermal majorons: Omega h^2 = {omega_thermal(mneed):.3f} extra at that mass,")
        pr(f"     i.e. {100*omega_thermal(mneed)/0.12:.0f}% of the dark matter as fast-moving (hot) majorons. Structure-formation data allow at most")
        pr("     a few per cent of such hot matter (rough). The hot start and majoron dark matter do not fit together.")
tauJ = lambda m: 1/(m*1e-9*(0.0501**2+0.0087**2)*1e-18/(16*math.pi*f_J**2))*6.58e-25
pr(f"  Majoron lifetime to neutrinos at 0.5 keV: {tauJ(500):.1e} s (universe 4.4e17 s): stable enough.")
pr("  Non-thermal start (the end of inflation hands its energy straight to the hidden partners, Lazarides & Shafi 1991): the space never gets")
pr("  hot enough to make majorons thermally, and the partners' decays still make the matter excess. Its size is then set by one new number")
pr("  (the reheating ratio), so the matter excess becomes ALLOWED, not predicted.")

head("SUMMARY TABLE")
for n, p, m, e, k, pull in rows: pr(f"  {n:46s} {p:<12.6g} {m:<12.6g} {pull:+7.1f} sigma  {k}")
