"""Why would the bag's inner space be octonionic? (14 Nov 2026)
Bag rules turned into algebra:
  Rule A (full glass, no gaps): two real relationships never combine into nothing  -> no zero divisors (ab=0 => a=0 or b=0)
  Rule B (conservation, water balloon): combining keeps size                       -> |ab| = |a||b|
Hurwitz (1898): the only real algebras with both are R, C, H, O (dims 1, 2, 4, 8).
Forces = symmetries of the inner space (the earlier rule, extended). Symmetry groups: R: none, C: a flip, H: SO(3), O: G2.
Tests: build all Cayley-Dickson algebras up to 16 (sedenions); check rules A and B; compute each symmetry algebra's size;
stabiliser of one direction in O; left/right action on H.
Credit: Hurwitz; Cayley; Dickson; Cartan (G2, 1914); Gunaydin & Gursey (1973); Dixon; Furey; Baez (2002); Furey & Hughes (2022)."""
import math, itertools, numpy as np
rng = np.random.default_rng(24)
def conj(x): return np.concatenate([x[:1], -x[1:]])
def mul(a, b):
    n = len(a)
    if n == 1: return a*b
    h = n//2; p, q, r, s = a[:h], a[h:], b[:h], b[h:]
    return np.concatenate([mul(p, r) - mul(conj(s), q), mul(s, p) + mul(q, conj(r))])
names = {1: "real numbers R", 2: "complex C", 4: "quaternions H", 8: "octonions O", 16: "sedenions S"}
def basis(n, i): e = np.zeros(n); e[i] = 1; return e

print("="*78); print("1. Which number systems obey the bag's two rules?"); print("="*78)
for n in (1, 2, 4, 8, 16):
    # Rule B: norm composition on random elements
    errB = max(abs(np.linalg.norm(mul(a, b)) - np.linalg.norm(a)*np.linalg.norm(b))
               for a, b in ((rng.normal(size=n), rng.normal(size=n)) for _ in range(2000)))
    # Rule A: search for zero divisors among (e_i +/- e_j)(e_k +/- e_l)
    zd = None
    if n >= 4:
        for i, j, k, l in itertools.product(range(1, n), repeat=4):
            if i < j and k < l:
                for s1, s2 in ((1, 1), (1, -1)):
                    a = basis(n, i) + s1*basis(n, j); b = basis(n, k) + s2*basis(n, l)
                    if np.abs(mul(a, b)).max() < 1e-12: zd = (i, j, s1, k, l, s2); break
            if zd: break
    tagA = "PASS (no zero divisors found)" if zd is None else f"FAILS: (e{zd[0]}{'+' if zd[2]>0 else '-'}e{zd[1]})(e{zd[3]}{'+' if zd[5]>0 else '-'}e{zd[4]}) = 0"
    tagB = "PASS" if errB < 1e-9 else f"FAILS (error {errB:.2f})"
    print(f"  dim {n:2d} {names[n]:16s} Rule A (no holes): {tagA:45s} Rule B (size kept): {tagB}")
print("  Hurwitz's theorem: exactly R, C, H, O pass, and nothing bigger can. The sedenions break both rules.")

print("\n"+"="*78); print("2. What forces can each one carry? (size of its symmetry algebra = number of carriers)"); print("="*78)
def derivations(n):
    E = [basis(n, i) for i in range(n)]
    T = np.array([[mul(E[i], E[j]) for j in range(n)] for i in range(n)])  # T[i,j,:] = e_i e_j
    rows = []
    for i in range(n):
        for j in range(n):
            for c in range(n):
                # D(e_i e_j)_c - (D e_i e_j)_c - (e_i D e_j)_c = 0 ; D unknown n x n, D[a,b] = (D e_b)_a
                row = np.zeros((n, n))
                for k in range(n): row[c, k] += T[i, j, k]
                for a in range(n):
                    row[a, i] -= T[a, j, c]
                    row[a, j] -= T[i, a, c]
                rows.append(row.ravel())
    M = np.array(rows); sv = np.linalg.svd(M, compute_uv=False)
    null = int((sv < 1e-9).sum()) + max(0, n*n - len(sv))
    return null, M
for n in (2, 4, 8, 16):
    d, _ = derivations(n)
    print(f"  {names[n]:16s} symmetry carriers: {d}")
print("  C: 0 continuous (only a flip). H: 3 = SO(3), the weak force's size. O: 14 = G2. S: still 14 (nothing new).")
print("  The strong force has 8 carriers that pull on each other: too big for H (3). Only O's symmetry can hold it.")

print("\n"+"="*78); print("3. Pick one direction in O (one 'perspective'): what is left?"); print("="*78)
n = 8; d, M = derivations(n)
U, S, Vt = np.linalg.svd(M); D = [Vt[-k-1].reshape(n, n) for k in range(d)]
for fix in range(1, 8):
    # derivations killing e_fix
    A = np.array([Dk[:, fix] for Dk in D]).T     # 8 x 14
    sv = np.linalg.svd(A, compute_uv=False); kept = d - int((sv > 1e-9).sum())
    print(f"  fix e{fix}: symmetries left = {kept}", end="")
    if fix == 7:
        # build the 8-dim stabiliser and get its rank via a generic element's centraliser
        _, s2, vt2 = np.linalg.svd(A); null = vt2[-kept:]
        G = [sum(c*Dk for c, Dk in zip(v, D)) for v in null]
        X = sum(rng.normal()*g for g in G)
        comm = np.array([(X@g - g@X).ravel() for g in G]).T
        rank = kept - np.linalg.matrix_rank(comm, tol=1e-8)
        # how it acts on the 6 directions orthogonal to 1 and e7: complex structure J = left mult by e7
        J = np.array([mul(basis(8, 7), basis(8, b)) for b in range(8)]).T
        ok = max(np.abs(g@J - J@g).max() for g in G)
        print(f"  (rank {rank}; commutes with 'multiply by e7' to {ok:.1e} -> acts on 3 complex directions) = SU(3)")
    else: print()
print("  Fixed direction + the real axis = 1 colourless complex slot; the other 6 = 3 colours.")

print("\n"+"="*78); print("4. Handedness inside H: left and right actions"); print("="*78)
q = rng.normal(size=4); q /= np.linalg.norm(q); x = rng.normal(size=4)
Lq = np.array([mul(q, basis(4, b)) for b in range(4)]).T; Rq = np.array([mul(basis(4, b), q) for b in range(4)]).T
print(f"  left- and right-multiplication by the same unit quaternion differ: |L - R| = {np.abs(Lq-Rq).max():.3f}")
print(f"  they commute with each other: |LR - RL| = {np.abs(Lq@Rq-Rq@Lq).max():.1e} -> two independent SU(2)s (left and right)")
print("  So H naturally carries SU(2)_L x SU(2)_R (Pati-Salam's two sides). WHICH side is light is not decided by H.")
