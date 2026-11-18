"""Neutron stars with nothing adjusted (14 Nov 2026).

Inward pull B comes from the proton (proton_sets_B_2026-11-14.py: B = 60.1, band 57-64 MeV/fm^3).
Outward push ('wall') comes from the known strong-force strength, not a knob:
  quark matter pressure  p(mu) = (3/(4 pi^2)) [a4 mu^4 - m_s^2 mu^2] - B,   a4 = 1 - 2 alpha_s/pi
  (Freedman & McLerran 1977; Farhi & Jaffe 1984; Fraga, Pisarski & Schaffner-Bielich 2001; Alford, Braby, Paris & Reddy 2005)
  alpha_s from two-loop running of the measured alpha_s(MZ) = 0.1180 down to scale Lbar = X mu, X = 2 to 4
  (the usual renormalisation-scale uncertainty, reported as a band, not tuned). m_s = 100 MeV; no pairing gap.
Tests: heaviest star >= 2.0 Msun (PSR J0740+6620, J0348+0432); Lambda(1.4) < ~580 (GW170817).
"""
import numpy as np, math, importlib.util, os
from scipy.interpolate import interp1d
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("tid", os.path.join(here, "neutron_star_tidal_fix_2026-11-13.py"))
tid = importlib.util.module_from_spec(spec); spec.loader.exec_module(tid)
hc = 197.3269

def alpha_s(Q):  # Q in GeV, two-loop, flavour thresholds at mb, mc
    def run(a, Q0, Q1, nf, n=400):
        b0 = (33-2*nf)/(12*math.pi); b1 = (153-19*nf)/(24*math.pi**2)
        t0, t1 = math.log(Q0**2), math.log(Q1**2); h = (t1-t0)/n
        f = lambda a: -(b0*a*a + b1*a**3)
        for _ in range(n):
            k1 = f(a); k2 = f(a+h*k1/2); k3 = f(a+h*k2/2); k4 = f(a+h*k3); a += h*(k1+2*k2+2*k3+k4)/6
        return a
    a = run(0.1180, 91.1876, 4.18, 5); a = run(a, 4.18, 1.27, 4); return run(a, 1.27, Q, 3)

def eos(B, X, ms=100.0, gap=0.0):
    mu = np.linspace(560.0/X, 1200, 3000)   # quark chemical potential, MeV
    a4 = np.array([1-2*alpha_s(X*m/1000)/math.pi for m in mu])
    p = (3/(4*math.pi**2))*(a4*mu**4 - (ms**2 - 4*gap**2)*mu**2)/hc**3 - B
    n = np.gradient(p, mu)               # quark density, fm^-3
    e = mu*n - p
    ok = p > 0
    i0 = np.argmax(ok)
    # surface: where p crosses 0
    es = np.interp(0.0, p[i0-1:i0+1], e[i0-1:i0+1])
    P, E = np.concatenate([[0.0], p[ok]]), np.concatenate([[es], e[ok]])
    return P, E, es, a4[i0]

def sequence(B, X, gap=0.0):
    P, E, es, a4s = eos(B, X, gap=gap)
    pf = interp1d(E, P, bounds_error=False, fill_value=(0.0, P[-1]))
    ef = interp1d(P, E, bounds_error=False, fill_value=(es, E[-1]))
    pf_ = lambda e: float(max(pf(e), 0.0))
    einv = lambda p: float(ef(p)) if p > 0 else es
    h = 0.5
    dpde = lambda e: max((pf_(e+h)-pf_(max(e-h, es)))/(e+h-max(e-h, es)), 1e-4)
    rows = np.array([tid.star(pf_, dpde, einv, ec, es_surf=es) for ec in np.geomspace(es*1.01, min(E[-1]*0.95, 3000), 70)])
    i = rows[:, 0].argmax(); r = rows[:i+1]
    L14 = float(np.exp(np.interp(1.4, r[:, 0], np.log(r[:, 3])))); R14 = np.interp(1.4, r[:, 0], r[:, 1])  # log(Lambda) interpolated: Lambda ~ C^-5
    return rows[i, 0], L14, R14, a4s, es

if __name__ == "__main__":
    print("alpha_s check: alpha_s(1.0 GeV) = %.3f, alpha_s(1.5 GeV) = %.3f" % (alpha_s(1.0), alpha_s(1.5)))
    print(" B (proton)  scale X   a4 at surface   Mmax (Msun)   R1.4 (km)   Lambda1.4   2.0 test   tidal test")
    for B in (57.1, 60.1, 63.7):
        for X in (2.0, 3.0, 4.0):
            M, L, R, a4, es = sequence(B, X)
            print(f"  {B:5.1f}      {X:3.1f}      {a4:.3f}          {M:.2f}         {R:5.2f}       {L:5.0f}      {'PASS' if M>=2.0 else 'fail'}       {'PASS' if L<580 else 'fail'}")

    print()
    print("What is missing: quark pairing (colour superconductivity; Alford, Rajagopal & Wilczek 1998), i.e. quarks")
    print("entangling in pairs at high density. Theory estimates the pairing gap at about 50-150 MeV; it is not measured.")
    print(" B      X    gap (MeV)   Mmax   R1.4   Lambda1.4   2.0 test   tidal test")
    for B in (57.1, 60.1, 63.7):
        for gap in (50.0, 100.0, 150.0):
            M, L, R, a4, es = sequence(B, 3.0, gap)
            print(f" {B:5.1f}  3.0   {gap:5.0f}      {M:.2f}   {R:5.2f}   {L:5.0f}      {'PASS' if M>=2.0 else 'fail'}       {'PASS' if L<580 else 'fail'}")
