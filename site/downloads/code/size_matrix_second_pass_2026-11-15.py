"""The second pass on the size matrix: (1) lay it out smallest to biggest as well as by family;
(2) line everything up by its equilibrium (balance) point with the proton, our home scale; (3) put the size and link matrices through the
eigen equations; (4) put the sizes through the cone (circulant eigen) equation the particles went through.
Credit: Koide 1981; Foot 1994; Brannen 2006; Harrison & Scott (circulants); Verhulst 1838 (logistic curve); PDG 2024."""
import math, os, numpy as np, openpyxl
from openpyxl.styles import Font
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
A = {"electron": 0.51099895, "muon": 105.6583755, "tau": 1776.86, "up": 2.16, "charm": 1273.0, "top": 172570.0,
     "down": 4.70, "strange": 93.5, "bottom": 4183.0, "nu1": 0.36e-9, "nu2": 8.66e-9, "nu3": 50.34e-9,
     "W": 80369.2, "Z": 91187.6, "Higgs": 125250.0, "proton": 938.272}
mp = A["proton"]
order = sorted(A, key=A.get); m = np.array([A[n] for n in order]); L = np.log(m); S = (L[:, None] - L[None, :])/2
head("1. Smallest to biggest: each particle's size from the proton (S = half log of mass/proton) and the step to the next")
for i, n in enumerate(order):
    step = (L[i] - L[i-1])/2 if i else float('nan')
    pr(f"  {n:9s} {m[i]:12.4g} MeV   S from proton {0.5*math.log(m[i]/mp):+8.3f}   step from previous {step:7.3f}")
steps = np.diff(L)/2
pr(f"  Biggest steps: {', '.join(f'{order[i]}->{order[i+1]} {steps[i]:.2f}' for i in np.argsort(steps)[::-1][:4])}")
pr("  The one big empty stretch is between the heaviest neutrino and the electron (S = 8.07, a factor 10 million in mass).")
pr("  Close neighbours from different families: muon and strange (0.06), W and Z (0.06; explained by the weak angle), proton, charm and tau")
pr("  (0.15, 0.17), Higgs and top (0.16). With 13 particles spread over S = -3.8 to +2.6 (average gap 0.53), about 1.4 gaps under 0.065")
pr("  are expected by chance; 2 are seen and one of them (W, Z) is already explained. Not a pattern by itself.")

head("2. Lined up by the equilibrium (balance) point with the proton")
pr("  The balance mass of any particle with the proton is m p/(m + p). As a fraction of the proton: phi = 1/(1 + exp(-2S)), with S the size from")
pr("  the proton. That is a logistic (S-shaped) curve: everything far lighter sits at its own mass, everything far heavier is capped at the proton,")
pr("  and the proton paired with itself sits exactly at one half.")
for n in order:
    s = 0.5*math.log(A[n]/mp); phi = 1/(1 + math.exp(-2*s)); M = 2*math.cosh(s)
    pr(f"  {n:9s} S {s:+8.3f}   balance with proton {phi*mp:12.5g} MeV ({phi:9.3e} of the proton)   link value {M:10.4g}")
pr("  What it shows: from the proton's home, the W, Z, Higgs and top all look almost the same (balance 0.988 to 0.993 of the proton), and the")
pr("  electron, neutrinos and light quarks are seen at their own masses. The proton's view compresses everything heavy into one cap. That is")
pr("  the scale-view idea of the record (K6) in exact form: each particle sees the others through this curve, centred on itself.")

head("3. The eigen equations applied to the size and link matrices themselves")
def link_eigs(names):
    x = np.sqrt(np.array([A[n] for n in names])); M = np.outer(x, 1/x) + np.outer(1/x, x)
    ev = np.sort(np.linalg.eigvalsh(M)); n = len(names); sm = sum(A[k] for k in names); si = sum(1/A[k] for k in names)
    return ev, n + math.sqrt(sm*si), n - math.sqrt(sm*si)
for label, names in (("9 charged matter particles", ["electron", "muon", "tau", "up", "charm", "top", "down", "strange", "bottom"]),
                     ("all 16", order)):
    ev, lp, lm = link_eigs(names)
    nz = ev[np.abs(ev) > 1e-6*abs(ev).max()]
    pr(f"  Link matrix, {label}: only {len(nz)} non-zero eigenvalues: {np.round(nz, 3).tolist()}")
    pr(f"     formula n +- sqrt(sum of masses x sum of 1/masses) = {lp:.3f}, {lm:.3f}")
ch = ["electron", "muon", "tau", "up", "charm", "top", "down", "strange", "bottom"]
sm = sum(A[k] for k in ch); si = sum(1/A[k] for k in ch)
pr(f"  For the charged particles: the top carries {100*A['top']/sm:.0f}% of the sum of masses; the electron {100/A['electron']/si:.0f}% of the sum of 1/masses")
pr(f"  (up quark {100/A['up']/si:.0f}%, down {100/A['down']/si:.0f}%). So the big eigenvalue, {9+math.sqrt(sm*si):.0f}, is mostly the electron-top link value ({2*math.cosh(0.5*math.log(A['top']/A['electron'])):.0f})")
pr(f"  stretched {math.sqrt(sm*si)/math.sqrt(A['top']/A['electron']):.2f} times by the light quarks. The matrix's eigen-structure is set mainly by its two ends.")
Lc = L - L.mean(); pr(f"  Size matrix (all 16): eigenvalues are 0 (14 times) and +-i x {0.5*math.sqrt(len(L)*np.sum(Lc**2)):.2f}: only the spread of the log masses survives.")
pr("  So the eigen equations confirm last night's muse in exact form: the matrix of relationships is governed by its lightest and heaviest member.")
pr("  They cannot produce the masses: both matrices are built from the masses and have only two directions (rank 2).")

head("4. The cone (circulant eigen) equation applied to the sizes instead of the square-root masses")
pr("  Particles: sqrt(mass) = a (1 + b cos(theta + 2 pi g/3)); for the leptons b = sqrt2 (the 45 degrees). Sizes: S_g = A + B cos(theta + 2 pi g/3), S from the proton.")
fam = {"leptons": ("electron", "muon", "tau"), "up type": ("up", "charm", "top"), "down type": ("down", "strange", "bottom"), "neutrinos": ("nu1", "nu2", "nu3")}
for f, ns in fam.items():
    s = np.array([0.5*math.log(A[n]/mp) for n in ns]); a0 = s.mean()
    c = np.sum((s - a0)*np.exp(-2j*math.pi*np.arange(3)/3)); B = 2*abs(c)/3; th = math.degrees(math.atan2(c.imag, c.real)) % 360
    r = np.sqrt([A[n] for n in ns]); r[0] *= (-1 if f == "neutrinos" else 1)
    a = r.mean(); cc = np.sum((r - a)*np.exp(-2j*math.pi*np.arange(3)/3)); b = 2*abs(cc)/3/a; th2 = math.degrees(math.atan2(cc.imag, cc.real)) % 360
    pr(f"  {f:10s} sizes: centre {a0:+7.3f}, swing {B:6.3f}, angle {th:7.2f} deg  |  sqrt-mass cone: swing/centre {b:.4f} (sqrt2 = 1.4142), angle {th2:7.2f} deg")
pr("  Every family fits (three numbers, three free values), so a fit alone says nothing. The sizes' angles differ family to family,")
pr("  so on sizes there is no shared cone. The square-root form stays the only one where a family (leptons, and neutrinos by construction) is exact.")

wb = openpyxl.Workbook(); ws = wb.active; ws.title = "Read_Me"
for r in ["Second pass on the size matrices (15 Nov 2026).", "Size_S_by_size: S = half log of mass ratio (row over column), particles smallest to biggest.",
          "Proton_Home: size from the proton, balance (equilibrium) mass with the proton, link value with the proton.",
          "The family-order sheets are in TRPG_Size_Matrices.xlsx."]: ws.append([r])
w = wb.create_sheet("Size_S_by_size"); w.append([""] + order)
for i, n in enumerate(order): w.append([n] + [round(float(v), 4) for v in S[i]])
for c in w[1]: c.font = Font(bold=True)
w = wb.create_sheet("Proton_Home"); w.append(["particle", "mass (MeV)", "S from proton", "balance mass with proton (MeV)", "balance / proton", "link value with proton"])
for n in order:
    s = 0.5*math.log(A[n]/mp); phi = 1/(1 + math.exp(-2*s)); w.append([n, A[n], round(s, 4), phi*mp, phi, 2*math.cosh(s)])
for c in w[1]: c.font = Font(bold=True)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TRPG_Size_Matrices_by_Size_and_Proton.xlsx"); wb.save(out); pr("\n  saved " + out)
