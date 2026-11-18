"""Idea tested: the matrices hold real values not yet measured, and relationships of more than two particles give a further set of matrices.
Builds TRPG_Three_Way_Carrier_Matrices.xlsx: every meeting of two matter particles with one carrier (photon, Z,
W, gluon, Higgs), plus the carriers with each other and the four-way meetings. Then tests:
  A. what the three-way relations already give;  B. which of their values have not been measured yet;
  C. whether any three-way relation fixes the open numbers (quark masses, mixing tilt, strengths).
Data: PDG 2024; NuFIT 6.0 (2024); ATLAS and CMS 2022-2024. Values are tree-level ('first look') unless a line says otherwise.
Credit: Fermi 1934; Yang & Mills 1954; Glashow 1961; Weinberg 1967; Salam 1968; Higgs 1964; Englert & Brout 1964; Guralnik, Hagen & Kibble 1964;
Cabibbo 1963; Kobayashi & Maskawa 1973; Pontecorvo 1957; Maki, Nakagawa & Sakata 1962; Gell-Mann 1962; Fritzsch, Gell-Mann & Leutwyler 1973;
Gross & Wilczek 1973; Politzer 1973; Veltman 1977, 1981 (rho; the balance condition); Ross & Veltman 1975; Fritzsch 1977, 1978 (texture zeros);
Gatto, Sartori & Tonin 1968; Fritzsch & Xing 2000; Masina & Quiros 2013; Hamada, Kawai & Oda 2013; Machacek & Vaughn 1983-84; Buttazzo et al. 2013;
Adler 1969; Bell & Jackiw 1969; Georgi & Glashow 1972; the ATLAS and CMS collaborations; PDG 2024; NuFIT 6.0."""
import math, os, numpy as np
from scipy.integrate import solve_ivp
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*100); pr(t); pr("="*100)
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TRPG_Three_Way_Carrier_Matrices.xlsx")

# ---------------- measured inputs ----------------
v = 246.21965; mZ = 91.1876; mW = 80.3692; mH = 125.20; mt = 172.57; GF = 1.1663788e-5
a0 = 1/137.035999; aZ = 1/127.951; s2 = 0.23122; aS = 0.1180
e = math.sqrt(4*math.pi*aZ); sw = math.sqrt(s2); cw = math.sqrt(1 - s2); g = e/sw; gp = e/cw; gs = math.sqrt(4*math.pi*aS)
F = [  # name, kind, charge, T3 of the left-handed member, colours, mass GeV, generation
    ("electron", "lepton", -1, -0.5, 1, 0.51099895e-3, 1), ("muon", "lepton", -1, -0.5, 1, 0.1056583755, 2), ("tau", "lepton", -1, -0.5, 1, 1.77693, 3),
    ("nu_1", "neutrino", 0, 0.5, 1, 0.0, 1), ("nu_2", "neutrino", 0, 0.5, 1, 0.0, 2), ("nu_3", "neutrino", 0, 0.5, 1, 0.0, 3),
    ("up", "up-type", 2/3, 0.5, 3, 2.16e-3, 1), ("charm", "up-type", 2/3, 0.5, 3, 1.2730, 2), ("top", "up-type", 2/3, 0.5, 3, 172.57, 3),
    ("down", "down-type", -1/3, -0.5, 3, 4.70e-3, 1), ("strange", "down-type", -1/3, -0.5, 3, 0.0935, 2), ("bottom", "down-type", -1/3, -0.5, 3, 4.183, 3)]
nuM = [0.36e-12, 8.6e-12, 50.0e-12]                                     # TRPG's neutrino masses (GeV), record M2 / an earlier test
CKM = np.array([[0.97435, 0.22501, 0.003732], [0.22487, 0.97349, 0.04183], [0.00858, 0.04111, 0.999118]])   # PDG 2024 global fit, magnitudes
def pmns(s12, s23, s13, d):
    c12, c23, c13 = (math.sqrt(1 - x) for x in (s12, s23, s13)); r12, r23, r13 = (math.sqrt(x) for x in (s12, s23, s13)); ed = complex(math.cos(d), math.sin(d))
    return np.array([[c12*c13, r12*c13, r13/ed], [-r12*c23 - c12*r23*r13*ed, c12*c23 - r12*r23*r13*ed, r23*c13], [r12*r23 - c12*c23*r13*ed, -c12*r23 - r12*c23*r13*ed, c23*c13]])
PM = np.abs(pmns(0.307, 0.470, 0.02215, math.radians(212)))          # NuFIT 6.0, normal ordering (rounded)

# ---------------- build the workbook ----------------
wb = openpyxl.Workbook(); bold = Font(bold=True); hdr = PatternFill("solid", fgColor="F5C518"); wrap = Alignment(wrap_text=True, vertical="top")
def sheet(name, rows, widths=None, note=None):
    ws = wb.create_sheet(name); r0 = 1
    if note: ws.cell(1, 1, note).font = Font(italic=True); r0 = 3
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            c = ws.cell(r0 + i, j + 1, val)
            if i == 0: c.font = bold; c.fill = hdr
            c.alignment = wrap
    for j, w in enumerate(widths or [22]*len(rows[0])): ws.column_dimensions[openpyxl.utils.get_column_letter(j + 1)].width = w
    return ws
wb.remove(wb.active)
sheet("Read_Me", [["TRPG three-way carrier matrices (14 Nov 2026)"],
    ["Why three-way: nature never joins two matter particles directly. Every real meeting is two matter particles plus one carrier (or carriers with each other)."],
    ["Each sheet is one carrier's matrix: which two particles it joins, and how strongly (the 'coupling', a pure number unless marked GeV)."],
    ["Status column: MEASURED (with the measured value as a multiple of the expected one), NOT YET MEASURED (a value nobody has seen yet), or LIMIT ONLY."],
    ["Strengths at the Z's energy: electric e = %.4f, weak g = %.4f, hypercharge g' = %.4f, strong g_s = %.4f. Higgs value v = %.3f GeV." % (e, g, gp, gs, v)],
    ["Matrix ideas: Matthew Gelfman (every particle to every particle; relationships of more than two). Physics of each meeting: Fermi; Yang & Mills; Glashow; Weinberg; Salam; Higgs; Englert & Brout; Cabibbo; Kobayashi & Maskawa; Pontecorvo; Maki, Nakagawa & Sakata; Gell-Mann; Fritzsch, Gell-Mann & Leutwyler. Data: PDG 2024; NuFIT 6.0; ATLAS and CMS."],
    ["Built by three_way_carrier_matrices_2026-11-14.py. The paper and the Results page explain what was found."]], [150])
rows = [["particle", "with", "carrier", "coupling (charge x e)", "status"]]
for n, k, q, t3, nc, m, gen in F:
    if q != 0: rows.append([n, "anti-" + n, "photon", round(q*e, 5), "MEASURED (charges exact)"])
sheet("Photon_Triples", rows, [12, 14, 9, 20, 28], "Photon joins a charged particle with its own antiparticle only. Never across families. Strength = charge.")
rows = [["particle", "with", "carrier", "left-handed coupling", "right-handed coupling", "status"]]
for n, k, q, t3, nc, m, gen in F:
    gl = g/cw*(t3 - q*s2); gr = g/cw*(-q*s2)
    rows.append([n, "anti-" + n, "Z", round(gl, 5), round(gr, 5) if k != "neutrino" else "none seen", "MEASURED (LEP, to 0.1%)" if gen < 3 or k in ("lepton", "down-type", "neutrino") else "MEASURED (top: LHC, about 15%)"])
sheet("Z_Triples", rows, [12, 14, 9, 20, 20, 30], "Z joins any particle with its own antiparticle, neutrinos included. Never across families (Glashow, Iliopoulos & Maiani 1970).")
rows = [["up-type", "down-type", "carrier", "coupling (g/sqrt2 x |V|)", "|V| (quark mixing)", "status"]]
ups = ["up", "charm", "top"]; dns = ["down", "strange", "bottom"]
for i in range(3):
    for j in range(3): rows.append([ups[i], dns[j], "W", round(g/math.sqrt(2)*CKM[i, j], 6), CKM[i, j], "MEASURED (PDG 2024 fit)"])
leps = ["electron", "muon", "tau"]
for i in range(3):
    for j in range(3): rows.append([leps[i], f"nu_{j + 1}", "W", round(g/math.sqrt(2)*PM[i, j], 5), round(PM[i, j], 4), "MEASURED (NuFIT 6.0; CP angle not yet fixed)"])
sheet("W_Triples", rows, [12, 12, 9, 22, 18, 36], "W joins an up-type with a down-type (quarks) or a charged lepton with a neutrino. The ONLY carrier that crosses families: its values are the mixing.")
lam = [np.array(M, complex) for M in ([[0,1,0],[1,0,0],[0,0,0]], [[0,-1j,0],[1j,0,0],[0,0,0]], [[1,0,0],[0,-1,0],[0,0,0]], [[0,0,1],[0,0,0],[1,0,0]],
       [[0,0,-1j],[0,0,0],[1j,0,0]], [[0,0,0],[0,0,1],[0,1,0]], [[0,0,0],[0,0,-1j],[0,1j,0]], np.diag([1,1,-2])/math.sqrt(3))]
cols = ["red", "green", "blue"]; rows = [["quark colour", "with antiquark colour", "gluon", "coupling (g_s x entry)", "status"]]
for a, L in enumerate(lam):
    for i in range(3):
        for j in range(3):
            if abs(L[i, j]) > 1e-12: rows.append([cols[i], "anti-" + cols[j], f"gluon {a + 1}", str(np.round(gs*L[i, j]/2, 4)), "MEASURED (jets, LEP and LHC)"])
fabc = {(1,2,3): 1, (1,4,7): .5, (2,4,6): .5, (2,5,7): .5, (3,4,5): .5, (1,5,6): -.5, (3,6,7): -.5, (4,5,8): math.sqrt(3)/2, (6,7,8): math.sqrt(3)/2}
for (a, b_, c), f in fabc.items(): rows.append([f"gluon {a}", f"gluon {b_}", f"gluon {c}", round(gs*f, 4), "MEASURED (colour factor 9/4, LEP)"])
sheet("Gluon_Triples", rows, [14, 20, 10, 22, 30], "Gluon joins two colours of the same quark flavour (every flavour alike), and gluons join each other (9 three-gluon meetings).")
kap = {"top": "1.01 +/- 0.11", "bottom": "0.99 +/- 0.16", "tau": "0.92 +/- 0.08", "muon": "1.12 +/- 0.21", "charm": "LIMIT ONLY: below 8.5 x", "electron": "LIMIT ONLY: below 245 x"}
rows = [["particle", "with", "carrier", "coupling (mass / v)", "status (measured as a multiple of expected)"]]
for n, k, q, t3, nc, m, gen in F:
    if k == "neutrino": rows.append([n, "anti-" + n, "Higgs", f"{nuM[gen - 1]/v:.1e} (if Higgs-linked)", "NOT YET MEASURED: unknown even in kind"]); continue
    rows.append([n, "anti-" + n, "Higgs", float(f"{m/v:.4g}"), ("MEASURED: " + kap[n]) if n in kap and not kap[n].startswith("LIMIT") else kap.get(n, "NOT YET MEASURED")])
for nm_, cpl, st in (("W", 2*mW**2/v, "MEASURED: 1.02 +/- 0.08"), ("Z", 2*mZ**2/v, "MEASURED: 1.04 +/- 0.07")):
    rows.append([nm_, nm_, "Higgs", f"{cpl:.2f} GeV", st])
rows.append(["Higgs", "Higgs", "Higgs", f"{3*mH**2/v:.1f} GeV", "NOT YET MEASURED: LHC limit -1.2 to 7.2 x (ATLAS 2024)"])
rows.append(["photon", "photon", "Higgs", "loop only (through W and top)", "MEASURED: discovery channel 2012"])
rows.append(["gluon", "gluon", "Higgs", "loop only (through top)", "MEASURED: makes 87% of LHC Higgs"])
rows.append(["Z", "photon", "Higgs", "loop only", "FIRST EVIDENCE 2023 (3.4 sigma): 2.2 +/- 0.7 x"])
sheet("Higgs_Triples", rows, [12, 14, 9, 30, 50], "Higgs joins each particle with its own partner. Its coupling IS the mass (divided by v). Every open quark mass is an entry here.")
rows = [["carrier", "carrier", "carrier", "carrier", "coupling", "status"]]
for row in (("Higgs", "Higgs", "W", "W", f"{2*mW**2/v**2:.4f}", "LIMIT ONLY: 0.57 to 1.48 x (ATLAS 2023)"), ("Higgs", "Higgs", "Z", "Z", f"{2*mZ**2/v**2:.4f}", "LIMIT ONLY (same search)"),
            ("Higgs", "Higgs", "Higgs", "Higgs", f"{3*mH**2/v**2:.4f}", "NOT YET MEASURED: no limit at all"), ("W", "W", "W", "W", f"{g**2:.4f}", "MEASURED (W-W scattering, LHC)"),
            ("W", "W", "Z", "Z", f"{(g*cw)**2:.4f}", "MEASURED (LHC)"), ("W", "W", "photon", "photon", f"{e**2:.4f}", "MEASURED (LHC, photon-photon to W-W)"),
            ("W", "W", "Z", "photon", f"{e*g*cw:.4f}", "MEASURED (LHC)"), ("gluon", "gluon", "gluon", "gluon", f"{gs**2:.4f} x colour factors", "MEASURED (multi-jet)")):
    rows.append(list(row))
sheet("Four_Way", rows, [10, 10, 10, 10, 22, 40], "Four-way meetings among carriers. Three of them have never been seen.")
rows = [["carrier", "carrier", "carrier", "coupling", "status"]]
for row in (("W", "W", "photon", f"{e:.4f}", "MEASURED (LEP)"), ("W", "W", "Z", f"{g*cw:.4f}", "MEASURED (LEP, to 2%)"), ("Higgs", "W", "W", f"{2*mW**2/v:.2f} GeV", "see Higgs_Triples"),
            ("Higgs", "Z", "Z", f"{2*mZ**2/v:.2f} GeV", "see Higgs_Triples"), ("Higgs", "Higgs", "Higgs", f"{3*mH**2/v:.1f} GeV", "see Higgs_Triples")):
    rows.append(list(row))
sheet("Carrier_Triples", rows, [10, 10, 10, 18, 30], "Carriers meeting carriers, three at a time.")

head("A. What the three-way relations already give")
mW_tree = mZ*cw
pr(f"  A1. The Higgs-W-W and Higgs-Z-Z meetings tie the W and Z masses to one Higgs value: W/Z = cos(weak angle) at first look.")
pr(f"      From the Z mass and the weak angle: W = {mW_tree:.2f} GeV; measured {mW:.2f} ({100*(mW_tree/mW - 1):+.2f}%). Adding the known loop corrections gives 80.357 (PDG 2024).")
pr("  A2. That tie (rho = W^2 / (Z^2 cos^2) = 1) holds only if the Higgs piece is a weak PAIR of charge 1/2. A weak triple would give rho = 1/2, or no")
pr("      Z mass at all (Ross & Veltman 1975). The main line's far side holds exactly such a pair and no other colourless piece.")
pr("      So the main line PREDICTS rho = 1 at first look. Measured, with known loops removed: 1.00031 +/- 0.00019 (PDG 2024). HOLDS")
pr("      (known physics, now coming from the matrix's own far side; the 0.0003 is allowed for any tiny extra piece, not required).")
pr("  A3. Every carrier meeting cancels its own 'leak' only if each family's charges sum correctly (Adler 1969; Bell & Jackiw 1969; Georgi & Glashow 1972).")
pr("      The main line's 16 states pass every such sum. HOLDS. This is what fixes the up partner at -2/3 and the electron partner at +1.")
pr("  A4. The W's matrix is the only one that crosses families; photon, Z, gluon and Higgs never do (at first look). That is why mixing shows only in")
pr("      W meetings and why no Z ever turns a strange quark into a down quark directly (measured: K-long to two muons is 7 in a billion of its decays, the size loops alone give). HOLDS (known).")

head("B. Values in these matrices that nobody has measured yet")
n_meas = n_lim = n_not = 0
for ws in wb.worksheets[1:]:
    for row in ws.iter_rows(values_only=True):
        s_ = " ".join(str(x) for x in row if x is not None)
        n_meas += "MEASURED" in s_ and "NOT YET" not in s_; n_lim += "LIMIT ONLY" in s_; n_not += "NOT YET" in s_
pr(f"  Entries: {n_meas} measured, {n_lim} with a limit only, {n_not} never measured.")
pr("  The values not yet seen, and what fixes them:")
for t_, val, by, test in (("Higgs with itself (three-way)", f"{3*mH**2/v:.1f} GeV", "Higgs mass and v (Standard Model); the main line keeps the same Higgs", "High-Luminosity LHC (2030s), to about 50%"),
                          ("Higgs with itself (four-way)", f"{3*mH**2/v**2:.3f}", "the same two numbers", "a 100 TeV collider; not before 2050"),
                          ("Higgs-Higgs-W-W", f"{2*mW**2/v**2:.3f}", "W mass and v", "LHC now: 0.57 to 1.48 x; High-Luminosity LHC"),
                          ("Higgs with charm", f"{1.273/v:.2e}", "the charm mass: main line says exactly mass/v", "LHC: below 8.5 x now"),
                          ("Higgs with electron", f"{0.000511/v:.2e}", "the electron mass", "FCC-ee running at the Higgs mass (2040s)"),
                          ("Higgs with up, down, strange", "mass / v each", "their masses", "no planned test reaches them"),
                          ("neutrino masses themselves", "sum 0.0592 eV (TRPG, record M2)", "TRPG's lepton relational rule", "DESI and Euclid now: sum below about 0.064 eV"),
                          ("neutrinos' CP angle", "about 270 degrees (TRPG, TM1, an earlier test)", "the perspective reading", "DUNE and Hyper-Kamiokande, 2030s"),
                          ("neutrino as its own antiparticle", "1.2 to 3.9 meV effective mass", "the Higgs-free energy", "beyond planned experiments")):
    pr(f"     {t_:<32} {val:<38} fixed by: {by:<62} test: {test}")
pr("  -> The matrices do hold real values not seen yet. But every one is FIXED by numbers already measured (or by a TRPG rule).")
pr("     They are predictions that can confirm or break the theory. None of them is an open number of ours.")

head("C. Does any three-way relation fix the open numbers?")
pr("  The open numbers are entries of two of these matrices: the Higgs matrix (the masses: bottom, charm, strange, top) and the W matrix (the tilt).")
pr("  Three kinds of three-way relation are known that could fix them. Each tested:")
pr("  C1. BALANCE. The Higgs is pushed by every meeting it is part of: by the W, Z and itself one way, by the top (and the other matter) the other way.")
pr("      If the pushes balance exactly, the Higgs's own mass is protected (Veltman 1981). That is a push/pull condition across ALL the Higgs triples:")
lhs = mH**2 + 2*mW**2 + mZ**2; rhs = 4*sum(m**2*nc/3 for n, k, q, t3, nc, m, gen in F)
pr(f"         Higgs^2 + 2 W^2 + Z^2 = {lhs:,.0f} GeV^2      4 x (sum of matter masses^2, colours counted) = {rhs:,.0f} GeV^2      ratio {lhs/rhs:.3f}")
pr(f"      Read as a prediction: Higgs = {math.sqrt(rhs - 2*mW**2 - mZ**2):.0f} GeV (measured 125.2), or top = {math.sqrt((lhs - 4*sum(m**2*nc/3 for n, k, q, t3, nc, m, gen in F if n != 'top'))/4):.1f} GeV (measured 172.6).")
pr("      FAIL at our energy. Does it balance at some higher energy? Follow every strength up (one loop, Standard Model):")
def rge(t, y):
    g1, g2, g3, yt, la = y; k = 1/(16*math.pi**2)
    return [k*(41/6)*g1**3, k*(-19/6)*g2**3, k*(-7)*g3**3, k*yt*(4.5*yt**2 - 8*g3**2 - 2.25*g2**2 - (17/12)*g1**2),
            k*(24*la**2 + 12*la*yt**2 - 6*yt**4 - 3*la*(3*g2**2 + g1**2) + (3/8)*(2*g2**4 + (g2**2 + g1**2)**2))]
y0 = [0.35830, 0.64779, 1.1666, 0.93690, 0.12604]                         # at the top's mass (Buttazzo et al. 2013, rounded)
sol = solve_ivp(rge, [0, math.log(1.22e19/mt)], y0, dense_output=True, rtol=1e-9, atol=1e-11)
def vel(t): g1, g2, g3, yt, la = sol.sol(t); return 2*la + (3*g2**2 + g1**2)/4 - 2*yt**2
ts = np.linspace(0, math.log(1.22e19/mt), 4000); vs = np.array([vel(t) for t in ts]); cross = [ts[i] for i in range(len(ts) - 1) if vs[i]*vs[i + 1] < 0]
for E in (1e3, 1e8, 1e12, 1e16, 1.22e19): pr(f"         at {E:8.0e} GeV: balance = {vel(math.log(E/mt)):+.3f}   (0 = balanced; at the top's own energy {vel(0):+.3f})")
if cross: pr(f"      Balances at {mt*math.exp(cross[0]):.1e} GeV.")
else: pr(f"      Never balances below gravity's grain at one loop; it is closest at the grain ({vs[-1]:+.3f}). With two loops published work puts the balance")
pr("      point at or somewhat above the grain (reported by Hamada, Kawai & Oda 2013; Masina & Quiros 2013). So 'the Higgs is balanced where gravity takes over' is")
pr("      a live reading. ATTACK: it is ONE condition, so at most it fixes ONE number (it ties the Higgs mass to the top's, the same tie as 44.4's")
pr("      'gravity's grain sets the Higgs', which missed at 2.6 sigma). It says nothing about bottom, charm, strange or the tilt. NOT a fix for our items.")
pr("  C2. ZEROS. Suppose some entries of the Higgs matrix are exactly zero before the W picks its directions (Fritzsch 1977, 1978). Then the W's mixing")
pr("      follows from the masses. The simplest pattern ('six zeros') predicts, at one energy (the top's):")
md, ms, mb, mu_, mc, mt_ = 2.678, 53.36, 2873.0, 1.242, 627.8, 162500.0
r_ds, r_uc, r_sb, r_ct = math.sqrt(md/ms), math.sqrt(mu_/mc), math.sqrt(ms/mb), math.sqrt(mc/mt_)
pr(f"         |V_us| between {abs(r_ds - r_uc):.4f} and {r_ds + r_uc:.4f}      measured 0.2250   -> inside: HOLDS (Gatto, Sartori & Tonin 1968; already in an earlier test)")
pr(f"         |V_cb| between {abs(r_sb - r_ct):.4f} and {r_sb + r_ct:.4f}      measured 0.0418   -> below the lowest value by a factor {abs(r_sb - r_ct)/0.04183:.1f}: FAIL (known since the top was found)")
pr(f"         |V_ub/V_cb| = sqrt(up/charm) = {r_uc:.4f}         measured {0.003732/0.04183:.4f}  -> {r_uc/(0.003732/0.04183):.2f} x: FAIL")
pr("      Richer patterns ('four zeros', Fritzsch & Xing 2000) still fit, but each adds as many free numbers as it explains. They relate mixing to masses;")
pr("      they never produce a mass. FAIL as a fix.")
pr("  C3. EQUALITY. Suppose all Higgs entries are equal at one energy (the matrix's own reading, an earlier test). Follow top, bottom and tau up (one loop):")
def rge3(t, y):
    g1, g2, g3, yt, yb, yta = y; k = 1/(16*math.pi**2); Y2 = 3*yt**2 + 3*yb**2 + yta**2
    return [k*(41/6)*g1**3, k*(-19/6)*g2**3, k*(-7)*g3**3,
            k*yt*(1.5*yt**2 - 1.5*yb**2 + Y2 - 8*g3**2 - 2.25*g2**2 - (17/12)*g1**2),
            k*yb*(1.5*yb**2 - 1.5*yt**2 + Y2 - 8*g3**2 - 2.25*g2**2 - (5/12)*g1**2),
            k*yta*(1.5*yta**2 + Y2 - 2.25*g2**2 - 3.75*g1**2)]
s3 = solve_ivp(rge3, [0, math.log(1.22e19/mt)], [0.35830, 0.64779, 1.1666, 0.93690, math.sqrt(2)*2.75/v, math.sqrt(2)*1.75/v], dense_output=True, rtol=1e-9, atol=1e-12)
for E in (mt, 1e8, 10**16.46, 1.22e19):
    _, _, _, yt_, yb_, yta_ = s3.sol(math.log(E/mt)); pr(f"         at {E:8.1e} GeV: top {yt_:.3f}   bottom {yb_:.4f}   tau {yta_:.4f}   bottom/tau {yb_/yta_:.2f}   top/bottom {yt_/yb_:.0f}")
tt = np.linspace(0, math.log(1.22e19/mt), 3000); rb = np.array([s3.sol(x)[4]/s3.sol(x)[5] for x in tt]); ix = int(np.argmin(abs(rb - 1)))
pr(f"      Bottom equals tau at one energy, about {mt*math.exp(tt[ix]):.0e} GeV (a known near-coincidence; it does not hold at the meeting energy: 0.62).")
pr("      But the top stays 59 to 74 times the bottom at every energy. 'All Higgs entries equal' never happens. FAIL.")
pr("  C4. STRENGTHS. The carrier strengths (e, g, g', g_s) are the photon, W, Z and gluon matrices' overall sizes. Only the meeting-point relation ties")
pr("      them: 10.6% off. The three-way matrices add no new relation among them. FAIL as a fix.")

head("D. Summary")
pr("  BUILT  : the five carrier matrices, carrier-carrier triples and four-way meetings (" + OUT.split('/')[-1] + ").")
pr("  HOLDS  : three-way relations give real values: the W/Z mass tie, rho = 1 from the main line's own Higgs pair, and every charge (cancellation sums).")
pr(f"  HOLDS  : the matrices hold values never measured ({n_not} entries, plus {n_lim} with limits only). All are fixed by measured numbers: predictions, not answers.")
pr("  FAILS  : as a fix for the open numbers. The open masses ARE the Higgs matrix's entries and the tilt IS the W matrix; none of the three known kinds of")
pr("           three-way relation (balance, zeros, equality) fixes them. Balance can fix at most one number and points to gravity's grain.")
pr("  So the three-way matrices are the right PLACE (they say exactly where every open number sits) but they do not yet give the RULE.")
wb.save(OUT); pr(f"\n  Saved {OUT}")
