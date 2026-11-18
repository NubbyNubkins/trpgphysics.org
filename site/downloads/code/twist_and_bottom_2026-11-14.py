"""2/9 twist and bottom mass: last honest tests (14 Nov 2026).
Twist: is delta = Q/3 (Koide value / number of generations) a rule shared with quarks? Is there a rule linking twist to the
'self vs relationship' probability (rule 7)?  Bottom: does rule 5 ('pairs share one relationship matrix') hold for quark pairs?
Credit: Koide; Brannen; Rivero & Gsponer; Georgi & Jarlskog; Pati & Salam; PDG 2024."""
import math, numpy as np
def lean_twist_Q(m):
    s = np.sqrt(np.array(m, float)); Q = sum(m)/s.sum()**2
    lean = math.degrees(math.acos(s.sum()/(math.sqrt(3)*np.linalg.norm(s))))
    z = sum(s[k]*np.exp(-2j*math.pi*k/3) for k in range(3)); d = math.atan2(z.imag, z.real) % (2*math.pi/3)
    return lean, d, Q
fams = {"leptons": [0.51099895, 105.6583755, 1776.93], "down (d,s,b)": [4.70, 93.5, 4183], "up (u,c,t)": [2.16, 1273, 172570]}
print("1. Twist rules")
for k, mm in fams.items():
    L, d, Q = lean_twist_Q(mm); p = math.cos(math.radians(L))**2
    print(f"  {k:14s} twist {d:.4f} | Q/3 = {Q/3:.4f} | (2/9)(2p) = {2/9*2*p:.4f} | (2/9)(1-Q)*3 = {2/9*3*(1-Q):.4f}")
print("  Q/3 matches leptons only (it IS 2/9 there by construction of Q = 2/3); no tested rule fits the quarks.")
print("\n2. Rule 5 on quark pairs: does the up member equal the down member (same Dirac matrix, no seesaw)?")
for (u, d, nm) in [(2.16, 4.70, "u/d"), (1273, 93.5, "c/s"), (172570, 4183, "t/b")]:
    print(f"  {nm}: ratio {u/d:8.2f} (rule 5 needs 1)")
print("  Ratios 0.46, 13.6, 41: not constant, so not even 'shared up to one factor'. Rule 5 FAILS for quarks.")
print("  For leptons it can't be seen directly, because the seesaw hides the neutrino's Dirac link. So rule 5 must be")
print("  restricted to leptons (colourless pairs), the same split as rule 7 (45 deg holds only without colour).")
print("\n3. Bottom = tau at the right-side scale (4x4 relation), 10^11.7 GeV")
for lg, r in [(11, 0.76), (16, 0.63)]: pass
r117 = 0.76 + (0.63-0.76)*(11.7-11)/(16-11)
print(f"  m_b/m_tau at 10^11.7 GeV ~ {r117:.2f} (interpolated from the earlier 0.76 at 10^11 and 0.63 at 10^16): {100*(1-r117):.0f}% off -> fails")
