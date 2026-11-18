"""Finishing the open items (14 Nov 2026).
1. Why one direction is picked in the octonions: does EVERY direction give the same SU(3)? (G2/SU(3) = S^6)
2. Why the right side is heavy: Pati-Salam left-right breaking scale vs the neutrino seesaw scale; one-loop unification
   (SM below M_I; SU(4)xSU(2)LxSU(2)R above) and a rough proton-lifetime check.
3. Quark masses: charm with all masses at one common energy (M_Z); the bottom input.
Credit: Cartan; Gunaydin & Gursey; Pati & Salam; Mohapatra & Pati; Senjanovic & Mohapatra; Minkowski; Gell-Mann, Ramond & Slansky;
Yanagida; Mohapatra & Senjanovic; Georgi, Quinn & Weinberg (1974, running); Chang, Mohapatra, Gipson, Marshak & Parida (1985);
Deshpande, Keith & Pal (1993); Koide; Huang & Zhou (2021, running masses); PDG 2024."""
import math, numpy as np
rng = np.random.default_rng(25)
def conj(x): return np.concatenate([x[:1], -x[1:]])
def mul(a, b):
    n = len(a)
    if n == 1: return a*b
    h = n//2; p, q, r, s = a[:h], a[h:], b[:h], b[h:]
    return np.concatenate([mul(p, r) - mul(conj(s), q), mul(s, p) + mul(q, conj(r))])
def e(i, n=8): v = np.zeros(n); v[i] = 1; return v
E = [e(i) for i in range(8)]
T = np.array([[mul(E[i], E[j]) for j in range(8)] for i in range(8)])
rows = []
for i in range(8):
    for j in range(8):
        for c in range(8):
            row = np.zeros((8, 8))
            for k in range(8): row[c, k] += T[i, j, k]
            for a in range(8): row[a, i] -= T[a, j, c]; row[a, j] -= T[i, a, c]
            rows.append(row.ravel())
M = np.array(rows); _, S, Vt = np.linalg.svd(M); D = [Vt[-k-1].reshape(8, 8) for k in range(14)]

print("="*78); print("1. Does every choice of direction leave the same strong force?"); print("="*78)
for trial in range(6):
    u = np.concatenate([[0], rng.normal(size=7)]); u /= np.linalg.norm(u)    # random unit imaginary octonion
    A = np.array([Dk@u for Dk in D]).T
    _, s2, vt2 = np.linalg.svd(A); kept = 14 - int((s2 > 1e-9).sum()); null = vt2[-kept:]
    G = [sum(c*Dk for c, Dk in zip(v, D)) for v in null]
    X = sum(rng.normal()*g for g in G)
    rank = kept - np.linalg.matrix_rank(np.array([(X@g - g@X).ravel() for g in G]).T, tol=1e-8)
    Ju = np.array([mul(u, E[b]) for b in range(8)]).T
    ok = max(np.abs(g@Ju - Ju@g).max() for g in G)
    print(f"  random direction {trial}: symmetries left {kept}, rank {rank}, acts on 3 complex directions (err {ok:.0e}) -> SU(3)")
# orbit: the 14 symmetries move a direction through 14 - 8 = 6 dimensions = the whole 6-sphere of unit imaginary octonions
print(f"  directions reachable by turning: 14 - 8 = 6 = every direction on the 6-sphere (G2/SU(3) = S^6).")
print("  -> Which direction is picked makes no difference; only THAT one is picked matters (like the Higgs).")

print("\n"+"="*78); print("2. Why is the right side heavy? Left-right breaking scale vs neutrino masses"); print("="*78)
a1, a2, a3 = 1/59.00, 1/29.59, 0.1180      # GUT-normalised alpha_1 at M_Z (PDG 2024)
MZ = 91.1876
bSM = np.array([41/10, -19/6, -7])
def solve(b4, b2L, b2R, label):
    # 1/a(mu) = 1/a(MZ) - b t, t = ln(mu/MZ)/(2 pi). Unknowns tI, tU.
    # At M_I: 1/a4 = 1/a3(I), 1/a2L = 1/a2(I), 1/a2R = (5/3)(1/a1(I) - (2/5)/a4(I))
    inv = np.array([1/a1, 1/a2, 1/a3])
    # write 1/aX(U) = c0 + c1 tI + c2 tU
    def at_I(i): return (inv[i], -bSM[i])
    c4 = (inv[2], -bSM[2] + b4, -b4)
    c2L = (inv[1], -bSM[1] + b2L, -b2L)
    c1I = at_I(0); c3I = at_I(2)
    r0 = (5/3)*(c1I[0] - 0.4*c3I[0]); r1 = (5/3)*(c1I[1] - 0.4*c3I[1])
    c2R = (r0, r1 + b2R, -b2R)
    Amat = np.array([[c4[1]-c2L[1], c4[2]-c2L[2]], [c2L[1]-c2R[1], c2L[2]-c2R[2]]])
    rhs = -np.array([c4[0]-c2L[0], c2L[0]-c2R[0]])
    tI, tU = np.linalg.solve(Amat, rhs)
    MI = MZ*math.exp(2*math.pi*tI); MU = MZ*math.exp(2*math.pi*tU)
    aU = 1/(c4[0] + c4[1]*tI + c4[2]*tU)
    tau = MU**4/(aU**2*0.938**5)*6.582e-25/3.156e7   # years, order of magnitude only
    print(f"  {label:44s} M_I = 10^{math.log10(MI):5.2f} GeV, M_U = 10^{math.log10(MU):5.2f} GeV, 1/a_U = {1/aU:4.1f}, proton ~10^{math.log10(tau):4.1f} yr")
    return MI, MU
# beta coefficients above M_I: gauge -11/3 C2, Weyl fermions +2/3 T, complex scalars +1/3 T; 3 generations of (4,2,1)+(4b,1,2)
f4 = f2 = 4.0
MIa, MUa = solve(-44/3 + f4 + 3, -22/3 + f2 + 1/3, -22/3 + f2 + 1/3 + 20/3, "minimal: bidoublet + right triplet (10,1,3)")
MIb, MUb = solve(-44/3 + f4 + 6, -22/3 + f2 + 1/3 + 20/3, -22/3 + f2 + 1/3 + 20/3, "left-right symmetric: + left triplet (10,3,1)")
y_top, m3 = 0.99, 0.0501e-9   # GeV
MR = (y_top*246.22/math.sqrt(2))**2/m3
print(f"  Seesaw: heaviest right-handed neutrino for our m3 = 0.0501 eV with a top-like link: M_R = {MR:.1e} GeV = 10^{math.log10(MR):.2f}")
print(f"  In left-right models M_R = f x v_R with f <~ 1, so v_R >~ M_R.")
print(f"  Minimal case: breaking scale 10^{math.log10(MIa):.1f} vs needed 10^{math.log10(MR):.1f}: {'consistent' if MIa > 0.3*MR else 'too low'}"
      f"; symmetric case 10^{math.log10(MIb):.1f}: {'consistent' if MIb > 0.3*MR else 'too low (needs f ~ ' + format(MR/MIb,'.0f') + ')'}")
print("  LHC: W_R > ~6 TeV; Super-K: proton > 2.4e34 yr (e+ pi0).")

print("\n"+"="*78); print("3. Quark masses: charm at one common energy"); print("="*78)
def koide_c(mt, mb):
    a, b = math.sqrt(mt), math.sqrt(mb); B = -4*(a+b); C = 3*(mt+mb) - 2*(a+b)**2
    return ((-B - math.sqrt(B*B-4*C))/2)**2, ((-B + math.sqrt(B*B-4*C))/2)**2
def Q(ms): return sum(ms)/sum(math.sqrt(x) for x in ms)**2
sets = {"own-scale (chain used before)": (172.57, 4.183, 1.273),
        "all at M_Z (Huang & Zhou 2021)": (168.26, 2.839, 0.620)}
for k, (mt, mb, mc) in sets.items():
    roots = [r for r in koide_c(mt, mb) if r < mb]
    pc = roots[0] if roots else float('nan')
    print(f"  {k:32s}: Q(c,b,t) = {Q((mc,mb,mt)):.4f}  charm from 45-deg rule = {pc:.3f} GeV vs {mc:.3f}  ({100*(pc-mc)/mc:+.1f}%)")
print("  Top at M_Z from 'full Higgs link' y_t = 1: v/sqrt2 = 174.1 GeV vs running 168.3 (+3.5%); pole 172.6 (+0.9%).")
print("\nBottom input: the 4x4 route (m_b = m_tau at the breaking scale) gives m_b/m_tau = 0.63-0.76 there: fails.")
print("  No relation in the matrices fixes the overall down-family scale; it stays the one quark input.")
