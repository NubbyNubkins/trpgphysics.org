"""The reverse-engineering of the sizes. Every particle pair gets its SIZE value
(S = half the log of the mass ratio, and the link value M = 2 cosh S) instead of the tracked values; the new matrix is searched for patterns,
each checked against chance (random mass sets) and against a second mass footing.
Footing A: the usual quoted masses (PDG 2024: pole for leptons and top, MS-bar for the others). Footing B: all at the Z's energy
(Xing, Zhang & Zhou 2008/2012). Neutrinos: TRPG's predicted masses. Credit: Froggatt & Nielsen 1979; Wolfenstein 1983;
Koide 1981; Xing, Zhang & Zhou 2008; PDG 2024."""
import math, os, itertools, numpy as np
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
A = {"electron": 0.51099895, "muon": 105.6583755, "tau": 1776.86, "up": 2.16, "charm": 1273.0, "top": 172570.0,
     "down": 4.70, "strange": 93.5, "bottom": 4183.0, "nu1": 0.36e-9, "nu2": 8.66e-9, "nu3": 50.34e-9,
     "W": 80369.2, "Z": 91187.6, "Higgs": 125250.0, "proton": 938.272}
B = {"electron": 0.486570, "muon": 102.718, "tau": 1746.24, "up": 1.38, "charm": 638.0, "top": 172100.0,
     "down": 2.82, "strange": 57.0, "bottom": 2860.0}
types = {"leptons": ("electron", "muon", "tau"), "up type": ("up", "charm", "top"), "down type": ("down", "strange", "bottom")}
def table(m): return np.array([[math.log(m[p]) for p in types[t]] for t in types])   # 3 types x 3 generations, natural-log masses
head("1. The size matrix (footing A). S = half log of mass ratio; link value M = 2 cosh S = (m1 + m2)/sqrt(m1 m2)")
names = list(A); L = np.array([math.log(A[n]) for n in names])
Smat = (L[:, None] - L[None, :])/2
pr("  Key fact: every size value is a DIFFERENCE of two log masses, S_ij = (L_i - L_j)/2. So the whole matrix holds exactly as many")
pr("  independent numbers as the particle list (one per particle, less one for the unit). Any pattern in it is a pattern in the log masses.")
pr("  " + " "*9 + "".join(f"{n[:7]:>8s}" for n in names[:12]))
for i, n in enumerate(names[:12]): pr(f"  {n[:8]:8s} " + "".join(f"{Smat[i, j]:8.2f}" for j in range(12)))

def stats(T):
    f2 = (T[:, 1] - T[:, 0])/(T[:, 2] - T[:, 0])               # where the middle generation sits between first and third (0..1)
    st1, st2 = T[:, 1] - T[:, 0], T[:, 2] - T[:, 1]
    resid_rank = T - T.mean(1, keepdims=True) - T.mean(0, keepdims=True) + T.mean()
    return f2, st1, st2, np.sqrt((resid_rank**2).sum()/4)
rng = np.random.default_rng(11)
def control(T, fn, N=20000):
    """random tables: each type keeps its lightest and heaviest mass (its span); its middle mass is placed at random (log-uniform) between."""
    out = []
    for _ in range(N):
        R = T.copy(); R[:, 1] = T[:, 0] + rng.uniform(0, 1, 3)*(T[:, 2] - T[:, 0]); out.append(fn(R))
    return np.array(out)
for foot, m in (("A (usual quoted masses)", A), ("B (all at the Z's energy)", B)):
    head(f"2. Patterns, footing {foot}")
    T = table(m); f2, st1, st2, rr = stats(T)
    pr(f"  log-mass span first->third: leptons {T[0,2]-T[0,0]:.3f}, up {T[1,2]-T[1,0]:.3f}, down {T[2,2]-T[2,0]:.3f}")
    pr(f"  P1 'same ratios in every family' (m = A_type x B_generation): rms miss in log {rr:.3f} (a factor {math.exp(rr):.2f} in mass). FAIL if >> 0.")
    pr(f"  P2 'same shape, stretched' (middle sits at the same fraction of the span): fractions {np.round(f2, 3)} (leptons, up, down)")
    spread = f2.max() - f2.min(); c = control(T, lambda R: (lambda f: f.max() - f.min())(stats(R)[0]))
    pr(f"      spread {spread:.3f}; random tables with a spread this small: {100*np.mean(c <= spread):.0f}%")
    r = f2[0]/f2[1], f2[1]/f2[2]
    pr(f"  P3 'fractions shrink by a fixed factor' (leptons -> up -> down): ratios {r[0]:.3f}, {r[1]:.3f}")
    c = control(T, lambda R: (lambda f: abs(math.log((f[0]/f[1])/(f[1]/f[2]))) if min(f) > 0 else 9)(np.sort(stats(R)[0])[::-1]))
    obs = abs(math.log(r[0]/r[1])); pr(f"      mismatch {obs:.3f}; random tables matching this well (fractions put in size order, the fair way): {100*np.mean(c <= obs):.0f}%")
    q = math.sqrt(1.5); pr(f"      one reading of the factor (found after looking, so a guess): sqrt(3/2) = {q:.4f} gives up {f2[0]/q:.3f}, down {f2[0]/q**2:.3f}")
    pr(f"  P4 'second step / first step' (spiral vs cone): leptons {st2[0]/st1[0]:.3f}, up {st2[1]/st1[1]:.3f}, down {st2[2]/st1[2]:.3f}")
    # P5 powers of one number: m/m_heaviest = eps^n, integers n
    best = None
    for eps in np.linspace(0.02, 0.5, 4801):
        n = (T - T[:, 2:3])/math.log(eps); n = n[:, :2]; dev = np.abs(n - np.round(n)).max()
        if np.round(n).min() >= 1 and (best is None or dev < best[0]): best = (dev, eps, np.round(n).astype(int))
    pr(f"  P5 'every mass = heaviest x (one number)^integer': best number {best[1]:.4f}, powers {best[2].tolist()}, worst off-integer {best[0]:.3f}")
    def p5(R):
        b = 9
        for eps in np.linspace(0.02, 0.5, 241):
            n = (R - R[:, 2:3])/math.log(eps); n = n[:, :2]
            if np.round(n).min() >= 1: b = min(b, np.abs(n - np.round(n)).max())
        return b
    c = control(T, p5, 1500); pr(f"      random tables doing at least this well (coarser search): {100*np.mean(c <= best[0] + 0.01):.0f}%")
    pr(f"  P6 partners in one generation (up vs down, log ratio): {np.round(T[1] - T[2], 3)}; leptons vs down: {np.round(T[0] - T[2], 3)}")
head("3. The neutrinos as a fourth family (TRPG's predicted masses)")
Ln = [math.log(A[k]) for k in ("nu1", "nu2", "nu3")]
pr(f"  middle fraction {(Ln[1]-Ln[0])/(Ln[2]-Ln[0]):.3f} vs charged leptons {stats(table(A))[0][0]:.3f}: the same because both follow the 45-degree rule")
pr("  (the neutrino masses were derived with it), so this is not an independent pattern.")

head("4. Write the size workbook")
import openpyxl
from openpyxl.styles import Font
wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Read_Me"
for r in ["TRPG size matrices (15 Nov 2026). The reverse-engineering of the sizes.",
          "S = half the natural log of the mass ratio (row over column). Link value M = 2 cosh S = (m1+m2)/sqrt(m1 m2).",
          "Every S is a difference of two log masses, so the matrix holds one number per particle (less the unit). Patterns in it are patterns in the log masses.",
          "Footing A: usual quoted masses (PDG 2024). Footing B: all at the Z's energy (Xing, Zhang & Zhou 2008). Neutrinos: TRPG's predicted masses.",
          "Pattern tests and their chance controls: the output of size_matrix_patterns_2026-11-15.py."]: ws.append([r])
def sheet(title, mat, nm, fmt):
    w = wb.create_sheet(title); w.append([""] + nm)
    for i, n in enumerate(nm): w.append([n] + [round(float(x), 4) for x in mat[i]])
    for c in w[1]: c.font = Font(bold=True)
sheet("Size_S_footingA", Smat, names, "S"); sheet("LinkValue_footingA", 2*np.cosh(Smat), names, "M")
nb = list(B); LB = np.array([math.log(B[n]) for n in nb]); SB = (LB[:, None] - LB[None, :])/2
sheet("Size_S_footingB", SB, nb, "S"); sheet("LinkValue_footingB", 2*np.cosh(SB), nb, "M")
w = wb.create_sheet("Family_Table")
w.append(["family", "log m1", "log m2", "log m3", "span", "middle fraction", "step1 (S)", "step2 (S)", "step2/step1", "footing"])
for foot, m in (("A", A), ("B", B)):
    T = table(m); f2, s1, s2_, _ = stats(T)
    for i, t in enumerate(types): w.append([t, *[round(x, 4) for x in T[i]], round(T[i,2]-T[i,0], 4), round(f2[i], 4), round(s1[i]/2, 4), round(s2_[i]/2, 4), round(s2_[i]/s1[i], 4), foot])
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TRPG_Size_Matrices.xlsx"); wb.save(out); pr("  saved " + out)
