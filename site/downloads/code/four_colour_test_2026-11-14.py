"""
Can the 3x3 colour matrix be a 4x4? (Pati & Salam 1974: lepton number as the 4th colour)
Tests: (1) does it keep what works, (2) does it fix the open problems, (3) what it predicts.
"""
import numpy as np
from fractions import Fraction as F

print("A. Plain 4 colours for QUARKS (N = 4) against data")
alpha, mpi, fpi = 1/137.036, 134.977, 92.07
for N in (3, 4):
    G = (N/3)**2*alpha**2*mpi**3/(64*np.pi**3*fpi**2)*1e6
    R = N*(4/9+1/9+1/9+4/9+1/9)*1.04
    print(f"   N={N}: pi0 width {G:5.2f} eV (meas 7.80+/-0.12), R {R:4.2f} (meas ~3.85), gluons {N*N-1}, C_A/C_F {2*N*N/(N*N-1):.3f} (meas ~2.2)")
print("   -> 4 quark colours FAILS. The working 4x4 is 3 quark colours + the lepton as the 4th.\n")

print("B. The Pati-Salam 4x4: SU(4) x SU(2)_L x SU(2)_R")
BL = [F(1,3)]*3 + [F(-1)]
print(f"   B-L inside the 4x4 = diag{tuple(str(x) for x in BL)}, trace = {sum(BL)}  (traceless => charges locked together)")
print("   SU(4) carriers: 15 = 8 gluons + 1 (B-L boson Z') + 3 + 3bar (6 leptoquark bosons, charge +-2/3)")
# charge from Q = T3L + T3R + (B-L)/2
for nm, t3l, t3r, bl in [("u_L",.5,0,1/3),("d_L",-.5,0,1/3),("nu_L",.5,0,-1),("e_L",-.5,0,-1),
                         ("u_R",0,.5,1/3),("d_R",0,-.5,1/3),("nu_R",0,.5,-1),("e_R",0,-.5,-1)]:
    print(f"   {nm:5s} Q = {t3l+t3r+bl/2:+.3f}")
print("   proton charge + electron charge = (2/3+2/3-1/3) + (-1) =", F(2,3)+F(2,3)-F(1,3)-1, "-> atoms are exactly neutral by the group itself (the 3x3 version needs the anomaly argument for this)")
# anomalies: SU(4)^3 for 16 = (4,2,1)+(4bar,1,2) written left-handed
print("   SU(4)^3 anomaly: 2 x A(4) + 2 x A(4bar) = 2 - 2 = 0 -> consistent with no extra conditions\n")

print("C. Quark masses: what the 4x4 links (mass relations at the Pati-Salam scale)")
# one-loop SM running of yb/ytau (Arason et al. 1992 beta functions), from m_t to M
def run(M, steps=4000):
    t0, t1 = np.log(172.6), np.log(M)
    g1, g2, g3 = 0.3583, 0.6477, 1.1645           # g' , g, g_s at m_t (MSbar)
    yt, yb, yta = 0.9369, 2.75*np.sqrt(2)/246.22, 1.7463*np.sqrt(2)/246.22
    h = (t1-t0)/steps; k = 1/(16*np.pi**2)
    for _ in range(steps):
        Y2 = 3*yt**2+3*yb**2+yta**2
        dg1 = k*41/6*g1**3; dg2 = -k*19/6*g2**3; dg3 = -k*7*g3**3
        dyt = k*yt*(1.5*yt**2-1.5*yb**2+Y2-(17/12*g1**2+9/4*g2**2+8*g3**2))
        dyb = k*yb*(1.5*yb**2-1.5*yt**2+Y2-(5/12*g1**2+9/4*g2**2+8*g3**2))
        dyta = k*yta*(1.5*yta**2+Y2-(15/4*g1**2+9/4*g2**2))
        g1+=h*dg1; g2+=h*dg2; g3+=h*dg3; yt+=h*dyt; yb+=h*dyb; yta+=h*dyta
    return yb/yta
for M in (1e11, 1e13, 1e14, 1e16):
    print(f"   m_b/m_tau at {M:.0e} GeV = {run(M):.3f}   (minimal 4x4 says 1.000)")
r = 105.6584/(9*0.51100)
print(f"   Georgi-Jarlskog (the 4x4's own -3 entry): m_s/m_d = m_mu/(9 m_e) = {r:.1f}; lattice ~ 20.0 -> {100*(r/20.0-1):+.0f}%")
print("   up and charm are tied to the (unknown) Dirac neutrino masses -> the 4x4 gives NO new handle on the up quark or charm\n")

print("D. Predictions")
mD, mnu = 172.6e9, 0.0502   # eV; heaviest Dirac ~ top (4x4 relation), m3 from atmospheric splitting
print(f"   right-handed neutrino (seesaw) scale ~ m_t^2/m_3 = {mD**2/mnu/1e9:.1e} GeV")
print("   leptoquark gauge bosons: K_L -> mu e < 4.7e-12 requires M > roughly 1000 TeV (out of LHC reach)")
print("   proton decay: Pati-Salam gauge bosons do NOT decay the proton (consistent with Super-K > 2.4e34 yr)")
print("   neutron <-> antineutron oscillation possible if B-L breaks by 2 (Mohapatra & Marshak 1980): ESS/NNBAR, DUNE")
print("   neutrinoless double-beta decay (Majorana neutrinos): LEGEND, nEXO, KamLAND-Zen")
