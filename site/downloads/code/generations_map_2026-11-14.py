"""Generations map and the 'three ways the bag can be occupied' test (14 Nov 2026).
an earlier test. Standard numbers: PDG 2024 masses; NuFIT 6.0 (2024) mixing; MEG II (2025); LEP Z width (Janot & Jadach 2020).
Tests:
  A. Lepton lean and twist (Koide 1981; Foot 1994; Brannen 2006).
  B. 'Same wave, three placements' = a circulant (Z3) relation matrix. Show its eigenvalues are exactly Brannen's form,
     and that lean 45 deg <=> self-link weight equals total cross-link weight.
  C. What three placements predict for mixing (trimaximal: Cabibbo 1978, Wolfenstein 1978; TM1/TM2: Lam; Albright & Rodejohann 2009).
  D. 'Three different shapes' reading (point, line, sheet defects) vs identical charges.
  E. Excited-states reading: does one level law fit?
  F. Octonion relation matrices: Jordan identity holds only up to 3x3 (Jordan, von Neumann & Wigner 1934).
  G. How many Standard Model numbers belong to the generations.
"""
import math, numpy as np
rng = np.random.default_rng(14112026)
print("="*78); print("A. Lepton lean and twist"); print("="*78)
me, mmu, mtau, dtau = 0.51099895, 105.6583755, 1776.93, 0.09
def lean_twist(m):
    s = np.sqrt(np.array(m, float))
    cos_t = s.sum()/(math.sqrt(3)*np.linalg.norm(s)); lean = math.degrees(math.acos(cos_t))
    mu = s.mean(); a = math.sqrt(2)*math.tan(math.radians(lean))  # amplitude in units of mu (sqrt2 at 45deg)
    # sqrt(m_k) = mu(1 + a cos(delta + 2 pi k/3)); solve delta from Fourier component
    z = sum(s[k]*np.exp(-2j*math.pi*k/3) for k in range(3))
    delta = (math.atan2(z.imag, z.real)) % (2*math.pi/3)
    return lean, delta, mu, a
for mt in (mtau-dtau, mtau, mtau+dtau):
    L, d, mu, a = lean_twist([me, mmu, mt])
    print(f"m_tau={mt:8.2f}: lean={L:.5f} deg, twist={d:.6f} rad (2/9={2/9:.6f}), off/diag weight={math.tan(math.radians(L))**2:.5f}")
Q = (me+mmu+mtau)/(math.sqrt(me)+math.sqrt(mmu)+math.sqrt(mtau))**2
print(f"Koide Q = {Q:.7f} (2/3 = {2/3:.7f})")

print("\n"+"="*78); print("B. Same wave in three placements = circulant (Z3) relation matrix"); print("="*78)
w = np.exp(2j*math.pi/3)
P = np.roll(np.eye(3), 1, axis=0)       # move one place round the circle
for trial in range(3):
    a = rng.uniform(0.5, 2); b = rng.uniform(0.1, 1)*np.exp(1j*rng.uniform(0, 2*math.pi))
    C = a*np.eye(3) + b*P + np.conj(b)*P.T
    ev = np.sort(np.linalg.eigvalsh(C))
    pred = np.sort([a + 2*abs(b)*math.cos(np.angle(b) + 2*math.pi*k/3) for k in range(3)])
    F = np.array([[w**(j*k) for k in range(3)] for j in range(3)])/math.sqrt(3)
    D = F.conj().T @ C @ F
    print(f"trial {trial}: eigenvalues {np.round(ev,6)}; Brannen form {np.round(pred,6)}; "
          f"max match error {abs(ev-pred).max():.1e}; Fourier diagonalises (off-diag {abs(D-np.diag(np.diag(D))).max():.1e})")
print("So: any 3x3 relation matrix that is unchanged by moving every member one place round (Z3)")
print("has eigenvalues a + 2|b| cos(arg b + 2 pi k/3): Brannen's formula, with twist = arg b, lean set by |b|/a.")
print("Lean 45 deg <=> 3a^2 (self weight) = 6|b|^2 (cross weight): tan^2(lean) = cross/self.")
for name, m in [("leptons", [me, mmu, mtau]), ("down quarks (d,s,b)", [4.70, 93.5, 4183]), ("up quarks (u,c,t)", [2.16, 1273, 172570])]:
    L, d, mu, a = lean_twist(m)
    print(f"  {name:22s} lean {L:6.3f} deg  cross/self weight {math.tan(math.radians(L))**2:.3f}  twist {d:.4f}")
print("Content check: 3 placements alone fit ANY three masses (3 numbers: scale, lean, twist).")
print("The physics is in the values: leptons have cross = self exactly (to 1e-5) and twist 2/9.")

print("\n"+"="*78); print("C. What 'three placements' predicts for mixing"); print("="*78)
s12, e12 = 0.307, 0.012; s13 = 0.02195; c13 = 1 - s13
print(f"NuFIT 6.0 (normal order): sin^2 th12 = {s12} +/- {e12}, sin^2 th13 = {s13}")
print(f"Pure trimaximal (charged leptons Z3, neutrinos aligned): sin^2 th13 = 1/3 -> excluded ({1/3:.3f} vs {s13})")
print(f"Tribimaximal (Harrison, Perkins & Scott 2002): sin^2 th13 = 0 -> excluded by Daya Bay/RENO 2012")
tm2 = 1/(3*c13); tm1 = 1 - 2/(3*c13)
print(f"TM2 (middle column all 1/3): sin^2 th12 = {tm2:.4f} -> {(tm2-s12)/e12:+.1f} sigma")
print(f"TM1 (first column 2/3,1/6,1/6): sin^2 th12 = {tm1:.4f} -> {(tm1-s12)/e12:+.1f} sigma (survives)")
print("Quarks: if up and down types are BOTH Z3-circulant, both are diagonalised by the same Fourier matrix,")
F = np.array([[w**(j*k) for k in range(3)] for j in range(3)])/math.sqrt(3)
print("  so CKM = F^dagger F =", np.round(abs(F.conj().T@F), 6).tolist(), "= identity (no mixing).")
print("  Measured |Vud|,|Vus|,|Vcb|,|Vub| = 0.9737, 0.2243, 0.0410, 0.0038: close to identity. Leading picture right;")
print("  the small mixing must come from the two families' placements differing slightly (not predicted here).")

print("\n"+"="*78); print("D. 'Three different shapes' (point, line, sheet) reading"); print("="*78)
print("In 3D, a localised disturbance can be a point (codim 3), a line (codim 2) or a sheet (codim 1): exactly 3 kinds.")
print("But their energies scale as size^0, size^1, size^2 and they carry DIFFERENT conserved charges and spins")
print("(Kibble 1976; Vilenkin & Shellard 1994). Generations carry IDENTICAL charges and spin 1/2 in every family:")
gens = {"e": (-1, -0.5, 0, 0.5), "mu": (-1, -0.5, 0, 0.5), "tau": (-1, -0.5, 0, 0.5)}
print("  (Q, T3, colour, spin) for e, mu, tau:", gens, "-> identical.")
print("Verdict: 'three different shapes' FAILS. Whatever the three ways are, they must be the same shape.")

print("\n"+"="*78); print("E. Excited-states reading: one level law?"); print("="*78)
for name, m in [("leptons", [me, mmu, mtau]), ("down", [4.70, 93.5, 4183]), ("up", [2.16, 1273, 172570])]:
    k12 = math.log(m[1]/m[0])/math.log(2); k23 = math.log(m[2]/m[1])/math.log(1.5)
    g12 = m[1]/m[0]; g23 = m[2]/m[1]
    print(f"  {name:8s} power law n^k: k(1->2)={k12:5.2f}, k(2->3)={k23:5.2f};  ratios {g12:8.1f}, {g23:8.1f}")
print("No single law within a family, and different between families. A tower also implies a 4th level:")
print("  Z width N_nu = 2.9963 +/- 0.0074 (LEP), and a 4th heavy chiral family is excluded by Higgs rates (LHC).")
print("  An excited muon would decay mu -> e gamma: MEG II (2025) BR < 1.5e-13. Composite leptons: scale > ~10-20 TeV (LEP).")
print("Verdict: a simple excitation tower is disfavoured (not strictly excluded: models can hide it).")

print("\n"+"="*78); print("F. Octonion-valued relation matrices: only up to 3x3 are consistent"); print("="*78)
def cd_mul(a, b):  # Cayley-Dickson product, arrays length 2^n
    n = len(a)
    if n == 1: return a*b
    h = n//2; p, q = a[:h], a[h:]; r, s = b[:h], b[h:]
    conj = lambda x: np.concatenate([x[:1], -x[1:]])
    return np.concatenate([cd_mul(p, r) - cd_mul(conj(s), q), cd_mul(s, p) + cd_mul(q, conj(r))])
def oconj(x): return np.concatenate([x[:1], -x[1:]])
def herm(n):
    M = np.zeros((n, n, 8))
    for i in range(n):
        M[i, i, 0] = rng.normal()
        for j in range(i+1, n):
            x = rng.normal(size=8); M[i, j] = x; M[j, i] = oconj(x)
    return M
def mm(A, B):
    n = A.shape[0]; C = np.zeros_like(A)
    for i in range(n):
        for j in range(n):
            for k in range(n): C[i, j] += cd_mul(A[i, k], B[k, j])
    return C
jp = lambda A, B: 0.5*(mm(A, B) + mm(B, A))
for dim, name in [(4, "quaternion"), (8, "octonion")]:
    for n in (2, 3, 4):
        errs = []
        for _ in range(6):
            if dim == 8: A, B = herm(n), herm(n)
            else:
                A, B = herm(n), herm(n); A[..., 4:] = 0; B[..., 4:] = 0
            A2 = jp(A, A)
            errs.append(np.abs(jp(jp(A, B), A2) - jp(A, jp(B, A2))).max())
        print(f"  {name:10s} {n}x{n}: Jordan identity error = {max(errs):.1e}  -> {'consistent' if max(errs) < 1e-9 else 'FAILS'}")
print("Theorem (Jordan, von Neumann & Wigner 1934): octonionic Hermitian matrices form a consistent algebra only up to 3x3")
print("(the 27-dimensional Albert algebra). Proposed link to 3 generations: Dubois-Violette 2016; Todorov & Dubois-Violette 2018;")
print("Boyle 2020. Link to generations is a published proposal, not proven.")

print("\n"+"="*78); print("G. Standard Model numbers that belong to the generations"); print("="*78)
items = [("charged fermion masses", 9), ("quark mixing angles", 3), ("quark CP phase", 1),
         ("neutrino masses", 3), ("lepton mixing angles", 3), ("lepton CP phase", 1), ("Majorana phases (if Majorana)", 2)]
tot = sum(n for _, n in items)
for nm, n in items: print(f"  {nm:32s} {n}")
other = 3 + 2 + 1   # 3 gauge couplings, Higgs mass and vev, theta_QCD
print(f"  generation-related: {tot} (or {tot-2} if neutrinos are Dirac);  everything else: {other} (3 couplings, 2 Higgs, theta)")
print(f"  share: {tot}/{tot+other} = {100*tot/(tot+other):.0f}%")
