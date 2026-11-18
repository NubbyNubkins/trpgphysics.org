"""Energy matrices built on the two kinds of force and the two perspectives (push/pull = strong, electric; change of perspective = weak, gravity).
For every pair of particles: how much energy the relationship holds under each of the four. Then test the frame on the matrices' shapes.
Credit: Bohr 1913; Fermi 1934; Sakharov & Zeldovich 1966; De Rujula, Georgi & Glashow 1975; Appelquist & Politzer 1975; Fadin & Khoze 1987 (toponium);
Gasiorowicz & Rosner 1981; Dirac 1937 / Eddington (large numbers); Stoney 1881; Cabibbo 1963; Kobayashi & Maskawa 1973; PDG 2024; CMS 2025 (top pair threshold)."""
import math, numpy as np
from scipy.optimize import least_squares
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
al = 1/137.035999; me = 0.51099895; mmu = 105.6583755; mtau = 1776.86; Mp = 938.27209; Mn = 939.56542; MPl = 1.220890e22   # MeV
v = 246219.65; mH = 125200.0; GF = 1.1663788e-11   # MeV^-2
P = {'e': (me, -1), 'mu': (mmu, -1), 'tau': (mtau, -1), 'p': (Mp, 1)}
mu_red = lambda a, b: a*b/(a + b)
link = lambda a, b: (a + b)/math.sqrt(a*b)            # the link value (every-particle-to-every-particle table)

head("1. PUSH/PULL #1, ELECTRIC. Energy of each pair bound by its own charges: (alpha q1 q2)^2 x (pair's balance mass) / 2   [Bohr 1913]")
pr("  The pair's balance mass (reduced mass) is EXACTLY sqrt(m1 m2) / (your link value): m1 m2/(m1+m2) = sqrt(m1 m2) x sqrt(m1 m2)/(m1+m2).")
pr("  So the every-particle table's link column already IS the binding matrix, divided into the pair's geometric-mean mass. Check, all pairs:")
worst = max(abs(mu_red(a, b) - math.sqrt(a*b)/link(a, b))/mu_red(a, b) for a in (me, mmu, mtau, Mp) for b in (me, mmu, mtau, Mp)); pr(f"     largest difference over 16 pairs: {worst:.1e} (identity)")
rows = [("hydrogen (p, e)", Mp, me, 13.598434e-6, "measured 13.598 eV"), ("positronium (e+, e)", me, me, 6.8029e-6, "measured 6.803 eV"),
        ("muonium (mu+, e)", mmu, me, 13.5390e-6, "theory/measured 13.54 eV"), ("muonic hydrogen (p, mu)", Mp, mmu, 2.5283e-3, "2.53 keV (from measured 2P-1S lines)"),
        ("true muonium (mu+, mu)", mmu, mmu, None, "NEVER OBSERVED (searches planned)"), ("tauonium (tau+, tau)", mtau, mtau, None, "never observed"),
        ("tau-hydrogen (p, tau)", Mp, mtau, None, "never observed")]
pr(f"  {'pair':26s} {'link value':>10s} {'binding (lowest level)':>24s}   comparison")
for n, a, b, meas, note in rows:
    B = al**2*mu_red(a, b)/2; s = f"{B*1e6:.4f} eV" if B < 1e-3 else f"{B*1e3:.3f} keV"
    pr(f"  {n:26s} {link(a, b):10.4f} {s:>24s}   {note}" + (f"  (difference {100*(B/meas - 1):+.3f}%, the rest is known small corrections)" if meas else ""))

head("2. PUSH/PULL #2, STRONG. Quark pairs: energy = masses + colour-magnetic link between every pair   [Sakharov & Zeldovich 1966; DGG 1975]")
pr("  Mesons (quark + antiquark) and baryons (three quarks). Each pair of quarks adds a(s1.s2)/(m1 m2): the spins' relationship divided by the masses.")
S0, S1 = -0.75, 0.25
mes = [("pi", 'q', 'q', S0, 138.0), ("K", 'q', 's', S0, 495.6), ("rho", 'q', 'q', S1, 775.3), ("omega", 'q', 'q', S1, 782.7), ("K*", 'q', 's', S1, 893.6),
       ("phi", 's', 's', S1, 1019.5), ("D", 'q', 'c', S0, 1867.2), ("D*", 'q', 'c', S1, 2008.6), ("Ds", 's', 'c', S0, 1968.4), ("Ds*", 's', 'c', S1, 2112.2),
       ("eta_c", 'c', 'c', S0, 2983.9), ("J/psi", 'c', 'c', S1, 3096.9), ("B", 'q', 'b', S0, 5279.5), ("B*", 'q', 'b', S1, 5324.7), ("Bs", 's', 'b', S0, 5366.9),
       ("Bs*", 's', 'b', S1, 5415.4), ("eta_b", 'b', 'b', S0, 9398.7), ("Upsilon", 'b', 'b', S1, 9460.3)]
# baryons: list of (pair, s.s) for the three pairs
bar = [("N", [('q', 'q', -0.25)]*3, 'qqq', 938.9), ("Delta", [('q', 'q', 0.25)]*3, 'qqq', 1232.0),
       ("Lambda", [('q', 'q', -0.75), ('q', 's', 0), ('q', 's', 0)], 'qqs', 1115.7), ("Sigma", [('q', 'q', 0.25), ('q', 's', -0.5), ('q', 's', -0.5)], 'qqs', 1193.2),
       ("Sigma*", [('q', 'q', 0.25), ('q', 's', 0.25), ('q', 's', 0.25)], 'qqs', 1384.6), ("Xi", [('s', 's', 0.25), ('q', 's', -0.5), ('q', 's', -0.5)], 'qss', 1318.3),
       ("Xi*", [('s', 's', 0.25), ('q', 's', 0.25), ('q', 's', 0.25)], 'qss', 1533.4), ("Omega", [('s', 's', 0.25)]*3, 'sss', 1672.5),
       ("Lambda_c", [('q', 'q', -0.75), ('q', 'c', 0), ('q', 'c', 0)], 'qqc', 2286.5), ("Sigma_c", [('q', 'q', 0.25), ('q', 'c', -0.5), ('q', 'c', -0.5)], 'qqc', 2453.5),
       ("Sigma_c*", [('q', 'q', 0.25), ('q', 'c', 0.25), ('q', 'c', 0.25)], 'qqc', 2518.1), ("Omega_c", [('s', 's', 0.25), ('s', 'c', -0.5), ('s', 'c', -0.5)], 'ssc', 2695.2),
       ("Lambda_b", [('q', 'q', -0.75), ('q', 'b', 0), ('q', 'b', 0)], 'qqb', 5619.6), ("Sigma_b", [('q', 'q', 0.25), ('q', 'b', -0.5), ('q', 'b', -0.5)], 'qqb', 5813.1)]
K = ['q', 's', 'c', 'b']
def mfit(p): m = dict(zip(K, p[:4])); return np.array([m[x] + m[y] + p[4]*1e7*s/(m[x]*m[y]) for _, x, y, s, _ in mes])
def bfit(p): m = dict(zip(K, p[:4])); return np.array([sum(m[c] for c in q) + p[4]*1e7*sum(s/(m[x]*m[y]) for x, y, s in pr_) for _, pr_, q, _ in bar])
for title, f, data in (("mesons", mfit, mes), ("baryons", bfit, bar)):
    meas = np.array([d[-1] for d in data]); r = least_squares(lambda p: (f(p) - meas)/meas, [330, 500, 1600, 4900, 1.0]); fit = f(r.x)
    pr(f"\n  {title}: 5 numbers fitted (four quark home masses + one link strength) against {len(data)} particles")
    pr(f"    home masses: up/down {r.x[0]:.1f}, strange {r.x[1]:.1f}, charm {r.x[2]:.1f}, bottom {r.x[3]:.1f} MeV")
    pr("    " + "  ".join(f"{d[0]} {fit[i]:.0f}/{meas[i]:.0f}" for i, d in enumerate(data)))
    pr(f"    typical miss {100*np.sqrt(np.mean(((fit - meas)/meas)**2)):.2f}%, worst {100*np.max(np.abs((fit - meas)/meas)):.1f}% ({data[int(np.argmax(np.abs((fit - meas)/meas)))][0]})")
    if title == "mesons": mq_mes = r.x[0]
pr(f"\n  TRPG reading: the up/down home mass from mesons, {mq_mes:.1f} MeV, against one third of the proton, {Mp/3:.2f} MeV (and the lepton cone's size 313.85):")
pr(f"    difference {100*(mq_mes/(Mp/3) - 1):+.1f}%. From baryons it is ~360 MeV, so this is a reading of one fit, not a law.")
pr("  The strong force's own 'hydrogen' (colour charge 4/3 between quark and antiquark), toponium, the heaviest pair:")
for a_s in (0.13, 0.15, 0.17):
    B = (4/3*a_s)**2*mu_red(172570, 172570)/2; pr(f"    alpha_s = {a_s}: binding {B/1e3:.2f} GeV  (top pair threshold excess seen by CMS and ATLAS 2025; theory ~2 GeV, Fadin & Khoze 1987)")

head("3. PERSPECTIVE #1, WEAK. The W matrix: each entry turns one particle's perspective into another's; energy = rate of the turn   [Fermi 1934]")
pr("  Rate of a turn from i to j: G_F^2 |V_ij|^2 m_i^5 /(192 pi^3) x (number of open ways). Every entry uses the same G_F (one strength).")
tmu = 2.1969811e-6; Btau = 0.1782
ttau = tmu*(mmu/mtau)**5*Btau; pr(f"  tau lifetime from the muon's (same strength, same shape; uses the measured share 17.82% for the electron route): {ttau*1e15:.1f} fs, measured 290.3 fs ({100*(ttau/290.3e-15 - 1):+.1f}%)")
mt, mW = 172570., 80369.2; x = (mW/mt)**2
Gt = GF*mt**3/(8*math.pi*math.sqrt(2))*(1 - x)**2*(1 + 2*x); Gtq = Gt*(1 - 0.09)
pr(f"  top's width (top turns to bottom, |V_tb| = 1): {Gt/1e3:.2f} GeV, with the strong-force correction {Gtq/1e3:.2f} GeV; measured 1.42 (+0.19 -0.15) GeV")
pr("  These are perspective turns in the plain sense: the particle stops being itself and becomes another. That is what the W entries do.")

head("4. PERSPECTIVE #2, GRAVITY. Energy of each pair held by gravity, against what its charges hold   [Dirac 1937; Stoney 1881]")
pr("  Gravity's entry for a pair is m1 m2 / M_Planck^2 (no charge, only mass). Ratio to the electric entry alpha q1 q2:")
for n, (a, qa), (b, qb) in (("e, p", P['e'], P['p']), ("e, e", P['e'], P['e']), ("p, p", P['p'], P['p']), ("mu, mu", P['mu'], P['mu'])):
    pr(f"    {n:7s} gravity / electric = {a*b/MPl**2/(al*abs(qa*qb)):.3e}")
pr(f"  Mass at which gravity equals electric for a charge-1 pair: sqrt(alpha) M_Planck = {math.sqrt(al)*MPl/1e3:.3e} GeV (Stoney's mass, 1881).")

head("5. THE SHAPES. Test of the frame: what kind of matrix is each relationship?")
q = np.array([2/3, -1/3, 0, -1]); m_el = np.array([2.16, 4.70, 1e-7, me])            # up, down, neutrino, electron (one family)
mats = {"electric (q1 q2)": np.outer(q, q), "gravity (m1 m2)": np.outer(m_el, m_el), "Higgs link (m1 m2 / v^2)": np.outer(m_el, m_el)/v**2}
for n, Mx in mats.items(): pr(f"  {n:26s}: rank {np.linalg.matrix_rank(Mx/np.abs(Mx).max(), tol=1e-9)} (a product: every pair = my number x your number; nobody changes identity)")
Vckm = np.array([[0.97435, 0.22501, 0.003732], [0.22487, 0.97349, 0.04183], [0.00858, 0.04111, 0.999118]])
pr(f"  W matrix (quarks)          : rank 3, rows each sum (squared) to {np.round((Vckm**2).sum(1), 4)}: a rotation (it TURNS one into another)")
pr("  gluon (colour)             : rotation in colour (turns red into blue) AND, summed over colours, a product (the colour-blind total that binds)")
pr("  Higgs vs gravity, for elementary particles: same matrix up to one number,")
pr(f"    Higgs entry / gravity entry = M_Planck^2/(4 pi v^2) = {MPl**2/(4*math.pi*v**2):.3e} for EVERY elementary pair")
pr(f"    (but the Higgs link reaches only {197.327/mH:.1e} fm, a thousandth of a proton; gravity reaches without limit)")
pr("  Results of the shape test:")
pr("   HOLDS  weak = perspective: the W matrix is a rotation, the only family-crossing one; its entries turn identity. Push/pull matrices never do.")
pr("   HOLDS  strong + electric = push/pull: product-shaped (colour-blind part for strong), they hold energy between partners who stay themselves.")
pr("   FAILS  gravity = perspective, in matrix form: gravity's matrix is a product (m1 m2), the push/pull shape, and changes no identity.")
pr("          It has the same shape as the Higgs link, 2e32 times weaker per pair at short range. Its perspective side (entanglement, geometry) is not visible")
pr("          in an energy matrix at all; it would have to show in which pairs are linked, not in how much energy they hold.")
pr("   HOLDS  the Z: its matrix is a product of weak charges, so the weak force has a push/pull half (Z) and a perspective half (W).")

head("6. Summary")
pr("  Energy matrices reproduce known values well (hydrogen, positronium, muonium exact to known corrections; 14 baryons to 0.9% and 18 mesons to 2.5%, 5 numbers each;")
pr("  tau lifetime and top width from one strength). They predict pairs never seen (true muonium 1.41 keV, tauonium 23.7 keV). All of this is")
pr("  known physics. What it says about TRPG: the link value is exactly the binding factor, and the four forces sort into products (push/pull) and rotations")
pr("  (perspective turns), with gravity on the product side. No open item moved: the matrices take the masses and strengths as inputs.")
