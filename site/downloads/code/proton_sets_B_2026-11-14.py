"""Set the containment depth B from the proton and the Delta (14 Nov 2026), with no adjustable numbers.

the containment balance: the bag radius R is where inward pull (B) balances outward quark push.
Model pieces (all from established physics, credited):
  - MIT bag: Chodos, Jaffe, Johnson, Thorn, Weisskopf (1974); DeGrand, Jaffe, Johnson, Kiskis (1975)
    E = (4pi/3) B R^3 + (3 x0 - Z0)/R + colour-magnetic term   (x0 = 2.0428 for massless quarks)
  - centre-of-mass correction M = sqrt(E^2 - <P^2>): Donoghue & Johnson (1980)
  - pion cloud ('leakage'): cloudy bag model, Theberge, Thomas & Miller (1980); Chew & Low static couplings
      Sigma_N = -(3 f^2/(4 pi^2 m^2)) [ I_N + (32/25) I_Delta ],  Sigma_Delta = -(3 f^2/(4 pi^2 m^2)) [ (8/25) PV I_N + I_Delta ]
      form factor u(kR) = 3 j1(kR)/(kR); f^2/4pi = 0.0778 (measured pion-nucleon coupling)
Fits: unknowns B and alpha_c to the measured N (938.9) and Delta (1232) masses; each radius from its own balance.
Z0 (the old free zero-point term) is set to 0 in the TRPG run: its job is done by the CM correction and the cloud.
"""
import numpy as np, math
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, fsolve

hc = 197.3269; x0 = 2.04279; MPI = 138.0; DM = 293.1
f2 = 0.0778*4*math.pi
MN_obs, MD_obs = 938.92, 1232.0
MAG = 0.17650*8*0.25    # massless colour-magnetic integral x pair sum (N: -, Delta: +)

def u(kR): return 1.0 if kR < 1e-6 else 3*(math.sin(kR)-kR*math.cos(kR))/kR**3
def I(R, dE, pv=False):
    Rk = R/hc
    if not pv:
        return quad(lambda k: k**4*u(k*Rk)**2/(math.sqrt(k*k+MPI*MPI)*(math.sqrt(k*k+MPI*MPI)+dE)), 0, 4000, limit=400)[0]
    # principal value at omega = DM (Delta -> N pi is open)
    k0 = math.sqrt(DM**2-MPI**2)
    # 1/(w - DM) = (w + DM)/(k^2 - k0^2); write integrand as h(k)/(k-k0)
    h = lambda k: k**4*u(k*Rk)**2/math.sqrt(k*k+MPI*MPI)*(math.sqrt(k*k+MPI*MPI)+DM)/(k+k0)
    return quad(h, 0, 4000, weight='cauchy', wvar=k0, limit=400)[0]
def sigma(R, which="N"):
    c = -3*f2/(4*math.pi**2*MPI**2)
    return c*(I(R, 0) + 32/25*I(R, DM))

def sigmaD(R):
    c = -3*f2/(4*math.pi**2*MPI**2)
    return c*(8/25*I(R, -DM, pv=True) + I(R, 0))

def mass(R, B, ac, Z0, which, cm=True, cloud=True):
    sgn = -1 if which == "N" else +1
    E = 4*math.pi/3*B*R**3 + hc*(3*x0 - Z0)/R + sgn*MAG*ac*hc/R
    M = math.sqrt(max(E*E - (3*(x0*hc/R)**2 if cm else 0), 1.0))
    if cloud:
        M += sigma(R, "N") if which == "N" else sigmaD(R)
    return M
def balanced(B, ac, Z0, which, cm=True, cloud=True):
    # Containment balance at the wall: quark push against B sets R (the cloud sits outside the wall, as leakage,
    # and is added afterwards; cloudy bag model practice).
    r = minimize_scalar(lambda R: mass(R, B, ac, Z0, which, cm=cm, cloud=False), bounds=(0.3, 2.5), method="bounded", options={"xatol":1e-5})
    return mass(r.x, B, ac, Z0, which, cm=cm, cloud=cloud), r.x
def fit(Z0, **kw):
    def eqs(v):
        B, ac = abs(v[0]), abs(v[1])
        return [balanced(B, ac, Z0, "N", **kw)[0]-MN_obs, balanced(B, ac, Z0, "D", **kw)[0]-MD_obs]
    B, ac = fsolve(eqs, [57.0, 1.5])
    B, ac = abs(B), abs(ac)
    mN, RN = balanced(B, ac, Z0, "N", **kw); mD, RD = balanced(B, ac, Z0, "D", **kw)
    return B, ac, RN, RD, mN, mD
def core_rms(R):
    j0 = lambda x: math.sin(x)/x
    j1 = lambda x: math.sin(x)/x**2 - math.cos(x)/x
    rho = lambda r: (j0(x0*r/R)**2 + j1(x0*r/R)**2)*r*r
    return math.sqrt(quad(lambda r: rho(r)*r*r, 1e-9, R)[0]/quad(rho, 1e-9, R)[0])

if __name__ == "__main__":
    print("1) Check against the 1975 fit (Z0 = 1.84, no CM correction, no cloud): expect B^1/4 near 145 MeV, B near 57")
    B, ac, RN, RD, mN, mD = fit(1.84, cm=False, cloud=False)
    print(f"   B = {B:.1f} MeV/fm^3 (B^1/4 = {(B*hc**3)**0.25:.0f} MeV), alpha_c = {ac:.2f}, R_N = {RN:.2f} fm")
    print("2) TRPG run: Z0 = 0 (no free zero-point term), CM correction on, pion cloud (leakage) on")
    B, ac, RN, RD, mN, mD = fit(0.0)
    print(f"   B = {B:.1f} MeV/fm^3 (B^1/4 = {(B*hc**3)**0.25:.0f} MeV), alpha_c = {ac:.2f}")
    print(f"   proton bag radius R = {RN:.3f} fm, quark-core rms = {core_rms(RN):.3f} fm, cloud energy = {sigma(RN,'N'):.0f} MeV")
    print(f"   check: M_N = {mN:.1f}, M_Delta = {mD:.1f}")
    print("3) Sensitivity: pion-nucleon coupling +/- 5%")
    for s in (0.95, 1.05):
        f2 = 0.0778*4*math.pi*s
        Bs = fit(0.0)[0]
        print(f"   coupling x{s}: B = {Bs:.1f} MeV/fm^3")
