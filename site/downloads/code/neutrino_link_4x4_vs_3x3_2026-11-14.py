"""Neutrino link: 4x4 (SO(10)/Pati-Salam) vs 3x3 (TRPG generation matrix). 14 Nov 2026, an earlier test.
4x4 route (minimal renormalisable SO(10): one 10 and one 126 Higgs; Babu & Mohapatra 1993; Bajc, Senjanovic & Vissani 2003):
   M_d = H + F,  M_e = H - 3F,  M_u = r(H + sF),  M_nuD = r(H - 3sF),  M_R = (v_R/v_d126) F     (H, F: 3x3)
   => F = (M_d - M_e)/4, H = (3M_d + M_e)/4: the right-handed neutrino masses are FIXED by down-quark minus lepton.
   Diagonal (no-mixing) approximation; masses at 2e16 GeV (SM running; Xing, Zhang & Zhou 2008), signs scanned.
3x3 route (TRPG): the neutrino links to the Higgs through the SAME 3x3 relation matrix as its doublet partner, the
   electron family (same lean 45 deg, same twist); seesaw with circulant (Z3) matrices.
Measured: m3/m2 = 5.76 (record M2; sqrt of the splitting ratio 5.8)."""
import math, itertools, numpy as np
u = [0.48e-3, 0.235, 74.0]; d = [1.14e-3, 0.0216, 0.99]; e = [0.4696e-3, 0.0991, 1.685]   # GeV at 2e16 GeV
m_nu = [0.00036e-9, 0.0087e-9, 0.0501e-9]   # GeV (record M2)
target32, target21 = m_nu[2]/m_nu[1], m_nu[1]/m_nu[0]
print("="*78); print("1. The 4x4 route (SO(10) with 10 + 126)"); print("="*78)
results = []
for sg in itertools.product([1, -1], repeat=6):   # signs of e2,e3 relative to d; signs of u1,u2,u3; sign of e1
    se = [sg[5], sg[0], sg[1]]; su = [sg[2], sg[3], sg[4]]
    F = [(d[i] - se[i]*e[i])/4 for i in range(3)]; H = [(3*d[i] + se[i]*e[i])/4 for i in range(3)]
    A = np.array([[H[1], F[1]], [H[2], F[2]]]); b = np.array([su[1]*u[1], su[2]*u[2]])
    try: r, rs = np.linalg.solve(A, b)
    except np.linalg.LinAlgError: continue
    nuD = [r*H[i] - 3*rs*F[i] for i in range(3)]
    u1pred = r*H[0] + rs*F[0]
    light = [nuD[i]**2/abs(F[i]) for i in range(3)]          # up to the common factor v_d/v_R
    p32 = light[2]/light[1]; p21 = light[1]/light[0]
    results.append((abs(math.log(p32/target32)), sg, r, rs/r if r else 0, nuD, F, p32, p21, u1pred))
results.sort(key=lambda x: x[0])
seen = set()
for res in results[:12]:
    key = round(res[6], 3)
    if key in seen: continue
    seen.add(key)
    _, sg, r, s, nuD, F, p32, p21, u1p = res
    print(f"  signs {sg}: r={r:7.2f} s={s:+7.3f} | Dirac nu links (GeV) {abs(nuD[0]):.2e} {abs(nuD[1]):.2e} {abs(nuD[2]):.2e}"
          f" | m3/m2 pred {p32:8.2f} (meas {target32:.2f}), m2/m1 pred {p21:9.1f} (rec {target21:.0f}), up pred {abs(u1p)*1e3:.2f} MeV (0.48)")
best = results[0]
_, sg, r, s, nuD, F, p32, p21, _ = best
MR3 = nuD[2]**2/m_nu[2]
print(f"\n  Best sign choice: heaviest Dirac link {abs(nuD[2]):.1f} GeV (top {u[2]} GeV) -> M_R3 = {MR3:.2e} GeV = 10^{math.log10(MR3):.1f}")
print(f"  v_R between M_R3 (126 coupling = 1) and M_R3 x 174/|F33| = 10^{math.log10(MR3*174/abs(F[2])):.1f} GeV")
print(f"  Allowed by the 4x4 running: 10^11.7 (minimal) to 10^13.7 (symmetric).")
spread = [x[6] for x in results]
print(f"  Across all sign choices m3/m2 ranges {min(spread):.2g} to {max(spread):.2g}; the measured 5.76 needs a specific choice.")

print("\n"+"="*78); print("2. The 3x3 route: neutrino linked like its partner, the electron family"); print("="*78)
mE = [0.51099895e-3, 0.1056584, 1.77693]   # GeV, low scale (leptons barely run)
MR = [mE[i]**2/m_nu[i] for i in range(3)]
for i, nm in enumerate(["1 (e)", "2 (mu)", "3 (tau)"]):
    print(f"  generation {nm:7s}: Dirac link = {mE[i]*1e3:9.3f} MeV -> hidden partner M_R = {MR[i]:.2e} GeV = 10^{math.log10(MR[i]):.2f}")
def Q(ms, sgn=(1, 1, 1)): return sum(ms)/sum(s*math.sqrt(m) for s, m in zip(sgn, ms))**2
print(f"  Koide Q of the hidden partners: {Q(MR):.4f} (all +); {Q(MR,(-1,1,1)):.4f} (lightest root -)  [2/3 = 0.6667]")
# Algebra: if sqrt(m_l) obeys Q=2/3 and sqrt(m_nu) obeys Q=2/3 with one sign flip, M_R = m_l^2/m_nu. Check self/cross of M_R
s = np.array([math.sqrt(x) for x in MR]); L = math.degrees(math.acos(s.sum()/(math.sqrt(3)*np.linalg.norm(s))))
print(f"  hidden partners' lean = {L:.2f} deg (leptons 45.00)")
print(f"  Right side must break at >= M_R3 = 10^{math.log10(MR[2]):.2f} GeV; 4x4 running allows 10^11.7: "
      f"{'CONSISTENT' if MR[2] < 10**11.7 else 'too high'} (126-type coupling needed: {MR[2]/10**11.74:.2f})")
print(f"  Leptogenesis check: lightest partner 10^{math.log10(MR[0]):.1f} GeV; simple thermal leptogenesis needs >~10^9 (Davidson & Ibarra 2002)"
      f" -> {'OK' if MR[0] > 1e9 else 'too light for the simple version (resonant versions still possible)'}")
print("  Note: the electron-family link uses the Higgs side; a neutrino uses the conjugate side (like up quarks). Using the")
print("  electron family's matrix is the TRPG 3x3 assumption (same relationship matrix for both members of a doublet).")
# what if the 3x3 link is the UP-quark family instead (conjugate side)?
mU = [2.16e-3, 1.273, 172.57]
MRu = [mU[i]**2/m_nu[i] for i in range(3)]
print(f"  Compare: up-family link (conjugate side, 4x4-like): M_R = 10^{math.log10(MRu[0]):.1f}, 10^{math.log10(MRu[1]):.1f}, 10^{math.log10(MRu[2]):.1f} GeV"
      f" -> needs the right side >= 10^{math.log10(MRu[2]):.1f}: too high for 10^11.7")
