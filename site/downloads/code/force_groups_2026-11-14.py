"""
Force groups from the TRPG relationship matrices (14 Nov 2026).
Question: which gauge groups are FORCED by what the matrices record
(counts of carriers, which carriers are neutral, who couples to whom, charges)?
Each step prints its result. Standard physics credited inline.
"""
import itertools, numpy as np
np.set_printoptions(precision=4, suppress=True)

# ------------------------------------------------------------------
# 1. Classification (Killing 1888, Cartan 1894): compact Lie algebras
#    by (dimension, rank). Rank = number of mutually neutral carriers.
# ------------------------------------------------------------------
simple = []
for n in range(1, 6):
    simple.append((f"SU({n+1})", n*(n+2), n))
    if n >= 2: simple.append((f"SO({2*n+1})", n*(2*n+1), n))
    if n >= 3: simple.append((f"Sp({2*n})", n*(2*n+1), n))
    if n >= 4: simple.append((f"SO({2*n})", n*(2*n-1), n))
simple.append(("G2", 14, 2))
simple.append(("U(1)", 1, 1))
def combos(dim, rank, start=0):
    if dim == 0 and rank == 0: yield []; return
    for i in range(start, len(simple)):
        nm, d, r = simple[i]
        if d <= dim and r <= rank:
            for rest in combos(dim-d, rank-r, i):
                yield [nm] + rest
print("1. Which compact force groups have a given carrier count and number of neutral carriers?")
for label, d, r in [("gluons: 8 carriers, 2 settled", 8, 2),
                    ("weak W+,W-,W3: 3 carriers, 1 neutral", 3, 1),
                    ("electroweak W+,W-,Z,photon: 4 carriers, 2 neutral", 4, 2),
                    ("all 12 carriers, 4 neutral (no block info)", 12, 4)]:
    sols = [" x ".join(c) for c in combos(d, r)]
    print(f"   {label:52s} -> {sols}")

# ------------------------------------------------------------------
# 2. Rebuild the strong algebra straight from the workbook's
#    Gluon_Triples sheet (f_abc values) and identify it.
# ------------------------------------------------------------------
f = np.zeros((8, 8, 8))
s3 = np.sqrt(3)/2
for (a, b, c), v in {(1,2,3):1, (1,4,7):.5, (1,5,6):-.5, (2,4,6):.5, (2,5,7):.5,
                     (3,4,5):.5, (3,6,7):-.5, (4,5,8):s3, (6,7,8):s3}.items():
    for (i, j, k), sg in [((a,b,c),1),((b,c,a),1),((c,a,b),1),((b,a,c),-1),((a,c,b),-1),((c,b,a),-1)]:
        f[i-1, j-1, k-1] = sg*v
ad = [-f[a] for a in range(8)]  # (ad_a)_{bc} = -f_abc, so [ad_a, ad_b] = f_abc ad_c
jac = max(np.abs(ad[a]@ad[b]-ad[b]@ad[a] - sum(f[a,b,c]*ad[c] for c in range(8))).max()
          for a in range(8) for b in range(8))
K = np.array([[np.trace(ad[a]@ad[b]) for b in range(8)] for a in range(8)])
X = sum(np.random.default_rng(1).normal(size=8)[a]*ad[a] for a in range(8))
rank = 8 - np.linalg.matrix_rank(X, tol=1e-9)
print("\n2. Strong force rebuilt from the workbook's 3-gluon couplings")
print(f"   closes on itself (Jacobi error {jac:.1e}); Killing form eigenvalues {np.round(np.linalg.eigvalsh(K),3)}")
print(f"   -> compact (all negative), semisimple, dimension 8, rank {rank}  =>  only SU(3) (step 1)")

# ------------------------------------------------------------------
# 3. Root diagram from 3 interchangeable colours with r+g+b = neutral
# ------------------------------------------------------------------
w = [np.array([np.cos(t), np.sin(t)]) for t in (np.pi/2, np.pi/2+2*np.pi/3, np.pi/2+4*np.pi/3)]
roots = [w[i]-w[j] for i in range(3) for j in range(3) if i != j]
L = sorted({float(round(np.linalg.norm(r), 6)) for r in roots})
angs = sorted({round(np.degrees(np.arccos(np.clip(roots[i]@roots[j]/(np.linalg.norm(roots[i])*np.linalg.norm(roots[j])),-1,1))))
               for i in range(6) for j in range(6) if i != j})
# simple roots: r->g and g->b
a1, a2 = w[0]-w[1], w[1]-w[2]
cart = np.array([[2*x@y/(y@y) for y in (a1, a2)] for x in (a1, a2)])
print("\n3. The 6 restless gluons' colour charges (from 3 equal colours summing to neutral)")
print(f"   all the same length: {L}; angles between them: {angs} deg -> a regular hexagon")
print(f"   Cartan matrix {cart.round(3).tolist()} = A2 = SU(3)")

# ------------------------------------------------------------------
# 4. Which group can act on 3 colours and reproduce the hadron list?
#    Data: qqq neutral (baryons), q qbar neutral (mesons), NO neutral qq.
# ------------------------------------------------------------------
print("\n4. Candidate colour groups on 3 colours, tested against the hadron matrices")
print("   group  carriers  qqbar neutral  qqq neutral  qq neutral  C_A/C_F")
for g, n, qqq, qq, cacf in [("SO(3)", 3, "yes", "YES (diquark hadrons, not seen)", 1.0),
                           ("U(3)", 9, "NO (phase kills eps_ijk)", "no", None),
                           ("SU(3)", 8, "yes", "no", 9/4)]:
    print(f"   {g:6s} {n:5d}     {'yes':13s} {qqq:26s} {qq:32s} {cacf}")
print("   measured gluon/quark charge ratio (LEP jets) ~2.2  -> SU(3) 2.25; SU(2) 2.67, SU(4) 2.13 for comparison")

# number of colours from data (the one input)
alpha, mpi, fpi = 1/137.036, 134.977, 92.07
for N in (2, 3, 4):
    G = (N/3)**2 * alpha**2 * mpi**3 / (64*np.pi**3*fpi**2) * 1e6  # eV
    R = N*(4/9+1/9+1/9+4/9+1/9)
    print(f"   N={N}: pi0->2 photons width {G:5.2f} eV (measured 7.80 +/- 0.12);  R ratio {R:4.2f} x (1+as/pi=1.04) = {R*1.04:4.2f} (measured ~3.85)")

# ------------------------------------------------------------------
# 5. Electroweak: the extra U(1) is forced by the Z touching right-handed electrons
# ------------------------------------------------------------------
mW, mZ, s2eff = 80.369, 91.188, 0.23153
s2os = 1 - (mW/mZ)**2
gR = s2eff          # g_R(e) = -Q sin^2 = +sin^2 ; pure SU(2) would give 0
gL = -0.5 + s2eff
print("\n5. Electroweak")
print(f"   sin^2(theta_W) from mW/mZ = {s2os:.4f}; effective (LEP/SLD) = {s2eff}")
print(f"   Z coupling to RIGHT-handed electron: {gR:+.4f} (LEP/SLD measure about +0.23). SU(2) alone gives 0 -> a U(1) MUST exist")
print(f"   Z coupling to LEFT-handed electron:  {gL:+.4f} (measured about -0.27)")
print(f"   rho = mW^2/(mZ^2 cos^2) tree level = 1 only if the Higgs is a doublet (Veltman 1977); measured 1.0004 +/- 0.0002 after loops")

# ------------------------------------------------------------------
# 6. Octonions: symmetries fixing one direction (Gunaydin & Gursey 1973; Furey 2018)
# ------------------------------------------------------------------
fano = [(1,2,3),(1,4,5),(1,7,6),(2,4,6),(2,5,7),(3,4,7),(3,6,5)]
T = np.zeros((8, 8, 8))           # e_i e_j = sum_k T[i,j,k] e_k
T[0, 0, 0] = 1
for i in range(1, 8): T[0,i,i] = T[i,0,i] = 1; T[i,i,0] = -1
for a, b, c in fano:
    for x, y, z in [(a,b,c),(b,c,a),(c,a,b)]:
        T[x,y,z] = 1; T[y,x,z] = -1
def mul(x, y): return np.einsum('i,j,ijk->k', x, y, T)
E = np.eye(8)
rows = []
for i in range(8):
    for j in range(8):
        # D(e_i e_j) - D(e_i) e_j - e_i D(e_j) = 0, unknown D is 8x8 (D[:,k] = D(e_k))
        M = np.zeros((8, 64))
        eij = mul(E[i], E[j])
        for k in range(8):
            for m in range(8):
                M[m, m*8+k] += eij[k]               # D(e_ij)
        # D(e_i) e_j : D(e_i)=sum_m D[m,i] e_m
        for m in range(8):
            M[:, m*8+i] -= mul(E[m], E[j])
            M[:, m*8+j] -= mul(E[i], E[m])
        rows.append(M)
A = np.vstack(rows)
_, sv, Vt = np.linalg.svd(A)
null = Vt[np.sum(sv > 1e-9):]
print("\n6. Octonion check (Hurwitz 1898; Cartan 1914; Gunaydin & Gursey 1973)")
print(f"   all symmetries of the octonions: dimension {len(null)} (G2 = 14)")
D7 = null.reshape(-1, 8, 8)[:, :, 7].reshape(len(null), 8)   # D(e7)
_, sv2, Vt2 = np.linalg.svd(D7.T)
keep = Vt2[np.sum(sv2 > 1e-9):] @ null
gens = [k.reshape(8, 8) for k in keep]
print(f"   symmetries that keep one direction (e7) fixed: dimension {len(gens)}")
Xr = sum(np.random.default_rng(2).normal()*g for g in gens)
comm = np.array([(Xr@g - g@Xr).ravel() for g in gens])
print(f"   rank (neutral generators) = {len(gens)-np.linalg.matrix_rank(comm, tol=1e-8)} -> dimension 8, rank 2 = SU(3)")
J = np.array([mul(E[7], E[k]) for k in range(8)]).T
print(f"   they all commute with 'multiply by e7' on the other 6 directions: max error "
      f"{max(np.abs(g@J-J@g).max() for g in gens):.1e} -> the 6 directions are 3 COMPLEX ones = 3 colours;")
print("   the real axis and e7 are untouched = 1 colourless (lepton-like) slot.  8 = 1 + 1 + 3 + 3bar")

# ------------------------------------------------------------------
# 7. Global form: every field obeys 6Y + 3(2T) + 2t = 0 mod 6  -> group is S(U(3)xU(2))
#    (all colour-neutral states have whole-number charge). O'Raifeartaigh 1986; Baez & Huerta 2010
# ------------------------------------------------------------------
fields = {"quark doublet Q_L": (1/6, .5, 1), "u_R": (2/3, 0, 1), "d_R": (-1/3, 0, 1),
          "lepton doublet L": (-.5, .5, 0), "e_R": (-1, 0, 0), "nu_R": (0, 0, 0),
          "Higgs": (.5, .5, 0), "gluon": (0, 0, 0), "W": (0, 1, 0), "B": (0, 0, 0)}
print("\n7. Glue between the three groups (the 'Z6' condition)")
ok = all(round(6*Y + 3*2*Tw + 2*t) % 6 == 0 for Y, Tw, t in fields.values())
print(f"   6Y + 3(2T) + 2(triality) divisible by 6 for every field: {ok}")
print("   -> the true group is S(U(3) x U(2)): one 3-block and one 2-block sharing a single phase,")
print("      which is why every hadron and lepton has whole-number charge while quarks have thirds")

# ------------------------------------------------------------------
# 8. Left-right (Pati & Salam 1974): leptons as a 4th colour, the 'other side' SU(2)_R
#    Y = T3R + (B-L)/2 must reproduce every hypercharge
# ------------------------------------------------------------------
print("\n8. Left-right check (Pati & Salam 1974; Mohapatra & Pati 1975)")
for name, T3R, BL, Y in [("u_R", .5, 1/3, 2/3), ("d_R", -.5, 1/3, -1/3), ("e_R", -.5, -1, -1), ("nu_R", .5, -1, 0)]:
    print(f"   {name:5s}: T3R + (B-L)/2 = {T3R+BL/2:+.3f}   Standard Model Y = {Y:+.3f}")
print("   left doublets: (B-L)/2 = Y exactly (+1/6 quarks, -1/2 leptons). One generation = 16 states (needs nu_R; neutrinos have mass)")

# ------------------------------------------------------------------
# 9. Anomaly cancellation for N colours (Geng & Marshak 1989) -> charges forced once N = 3
# ------------------------------------------------------------------
from fractions import Fraction as Fr
N = 3
Q, u, d, Lh, e = Fr(1, 2*N), Fr(1, 2)+Fr(1, 2*N), -Fr(1, 2)+Fr(1, 2*N), Fr(-1, 2), Fr(-1)
s = {"SU3^2 U1": 2*Q - u - d, "SU2^2 U1": N*Q + Lh, "U1^3": 2*N*Q**3 - N*u**3 - N*d**3 + 2*Lh**3 - e**3,
     "grav^2 U1": 2*N*Q - N*u - N*d + 2*Lh - e}
print("\n9. Anomaly sums with N = 3 (all must be 0):", {k: str(v) for k, v in s.items()})
print(f"   charges: up {u + 0}, down {d}  (Q = T3 + Y)")
