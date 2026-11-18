"""The master mass tables. (1) Do the leftover numbers in the matrices (the balance points of
every pair) match real masses more often than chance? (2) Each mass against its equilibrium points. (3) Master mass table built from three
layers: the Higgs field (elementary masses) + the equilibrium (the strong force's containment balance, which turns
quarks into protons and other hadrons) + the energy of the nucleus's interactions (binding energy); tested against measured masses.
Credit: Weizsacker 1935; Bethe & Bacher 1936 (liquid drop); De Rujula, Georgi & Glashow 1975 (hadron masses); Sakharov & Zeldovich 1966;
Ji 1995 and Yang et al. (chiQCD) 2018 (how the proton's mass splits); Wang et al., AME2020 (nuclear masses, via the periodictable package);
PDG 2024."""
import math, sys, itertools, numpy as np
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
P = {"electron": 0.51099895, "muon": 105.6583755, "tau": 1776.86, "up": 2.16, "down": 4.70, "strange": 93.5, "charm": 1273.0, "bottom": 4183.0,
     "top": 172570.0, "W": 80369.2, "Z": 91187.6, "Higgs": 125250.0, "proton": 938.272, "neutron": 939.565, "pion": 139.570, "kaon": 493.677}
head("1. Leftover numbers: does the balance point of a pair land on another particle's mass more often than chance?")
names = list(P); m = np.array([P[n] for n in names])
def matches(m, tol=0.02):
    hits = []
    for i, j in itertools.combinations(range(len(m)), 2):
        mu = m[i]*m[j]/(m[i] + m[j])
        for k in range(len(m)):
            if k in (i, j): continue
            if abs(mu/m[k] - 1) < tol and abs(mu/min(m[i], m[j]) - 1) > tol: hits.append((i, j, k, mu))   # skip 'balance = the lighter one' cases
    return hits
h = matches(m)
for i, j, k, mu in h: pr(f"  balance({names[i]}, {names[j]}) = {mu:10.4g} MeV  ~  {names[k]} {m[k]:.4g} MeV ({100*(mu/m[k]-1):+.2f}%)")
rng = np.random.default_rng(5); lo, hi = np.log(m.min()), np.log(m.max()); cnt = []
for _ in range(3000): cnt.append(len(matches(np.exp(rng.uniform(lo, hi, len(m))))))
cnt = np.array(cnt); pr(f"  {len(h)} matches within 2% (the trivial 'balance = lighter partner' cases removed). Random mass lists of the same size and range:")
pr(f"  average {cnt.mean():.1f} matches; {100*np.mean(cnt >= len(h)):.0f}% have at least {len(h)}.")

head("2. Each mass against its equilibrium points")
v = 246.21965; hb = v/math.sqrt(2)
pr(f"  Higgs balance point v/sqrt2 = {hb:.2f} GeV. Each elementary mass as a fraction of it (the Higgs link):")
for n in ("electron", "up", "down", "strange", "muon", "charm", "tau", "bottom", "top"):
    pr(f"   {n:9s} {P[n]/1000/hb:10.3e}")
pr("  W, Z, Higgs: g v/2, sqrt(g^2+g'^2) v/2, sqrt(2 lambda) v: their balance comes from the force strengths and the Higgs self-coupling.")
pr("  With the proton (our home), every particle's balance mass is m p/(m + p) = p/(1 + exp(-2S)).")

head("3a. Layer 2, the equilibrium: quark home masses that make the hadrons (strong containment balance + spin interaction)")
B = {"N": 938.92, "Delta": 1232.0, "Lambda": 1115.683, "Sigma": 1193.15, "Xi": 1318.28, "Sigma*": 1384.6, "Xi*": 1533.4, "Omega": 1672.45,
     "Lambda_c": 2286.46, "Lambda_b": 5619.60}
def baryons(p):
    u, s, c, b, A = p; A *= 1e7
    return {"N": 3*u - 3*A/u**2, "Delta": 3*u + 3*A/u**2, "Lambda": 2*u + s - 3*A/u**2, "Sigma": 2*u + s + A*(1/u**2 - 4/(u*s)),
            "Xi": u + 2*s + A*(1/s**2 - 4/(u*s)), "Sigma*": 2*u + s + A*(1/u**2 + 2/(u*s)), "Xi*": u + 2*s + A*(1/s**2 + 2/(u*s)),
            "Omega": 3*s + 3*A/s**2, "Lambda_c": 2*u + c - 3*A/u**2, "Lambda_b": 2*u + b - 3*A/u**2}
from scipy.optimize import least_squares
fit = least_squares(lambda p: [(baryons(p)[k] - B[k])/B[k] for k in B], [360, 540, 1700, 5040, 0.5])
pb = fit.x; pred = baryons(pb)
for k in B: pr(f"   {k:9s} measured {B[k]:8.2f}  from home masses {pred[k]:8.2f}  ({100*(pred[k]/B[k]-1):+.2f}%)")
light = [k for k in B if k not in ("Lambda_c", "Lambda_b")]
rms = math.sqrt(np.mean([(pred[k]/B[k]-1)**2 for k in light]))
pr(f"  The charm and bottom home masses are each set by one baryon, so the real test is the 8 light baryons from 3 numbers: rms {100*rms:.2f}%.")
cur = {"up/down": (P["up"] + P["down"])/2, "strange": P["strange"], "charm": P["charm"], "bottom": P["bottom"]}
pr("  Home mass = Higgs part (current mass) + equilibrium part (what the strong containment adds):")
for (n, cm), hm in zip(cur.items(), pb[:4]): pr(f"   {n:8s} home {hm:7.1f} MeV = Higgs {cm:7.1f} + equilibrium {hm - cm:6.1f} MeV")
pr("  The equilibrium part is a few hundred MeV for every quark: the same containment energy, roughly, whatever the quark (charm and bottom")
pr("  depend on which mass definition is used: with the bottom's pole mass, 4.78 GeV, its equilibrium part is about 260 MeV).")
pr("  The proton itself, split by lattice calculations (Yang et al. 2018): quark masses (the Higgs layer) about 9%, quark motion 32%, gluon")
pr("  field 37%, the scale (trace) part 23%. So about 91% of the proton is equilibrium energy, as the record's containment principle says.")

head("3b. Layer 3, the nucleus's interaction energy: binding energy from the balance (liquid drop) formula, against AME2020")
import periodictable as pt
mu_amu = 931.49410242; me = P["electron"]; mH = 938.272088 + me - 13.6e-6; mn = 939.565421
rows = []
for el in pt.elements:
    if el.number < 1: continue
    for iso in el:
        A_ = iso.isotope; Z = el.number; N = A_ - Z
        if A_ < 12 or iso.mass is None or iso.mass <= 0: continue
        M = iso.mass*mu_amu; Bexp = Z*mH + N*mn - M
        rows.append((Z, N, A_, M, Bexp, f"{el.symbol}-{A_}"))
Z, N, Aa, M, Bx = (np.array([r[i] for r in rows], float) for i in range(5)); lab = [r[5] for r in rows]
delta = np.where((Z % 2 == 0) & (N % 2 == 0), 1, np.where((Z % 2 == 1) & (N % 2 == 1), -1, 0))/np.sqrt(Aa)
X = np.column_stack([Aa, -Aa**(2/3), -Z*(Z - 1)/Aa**(1/3), -(N - Z)**2/Aa, delta])
coef, *_ = np.linalg.lstsq(X, Bx, rcond=None); Bp = X @ coef
pr(f"  {len(rows)} nuclei with A >= 12. Balance formula: volume {coef[0]:.2f}, surface {coef[1]:.2f}, electric {coef[2]:.3f}, asymmetry {coef[3]:.2f}, pairing {coef[4]:.1f} MeV")
res = Bp - Bx; pr(f"  binding energy rms miss {np.sqrt(np.mean(res**2)):.2f} MeV (binding energies run up to {Bx.max():.0f} MeV)")
pr("  Master mass of each atom = Z x hydrogen atom + N x neutron - binding:")
Mp = Z*mH + N*mn - Bp; rel = (Mp - M)/M
pr(f"   without the binding layer: worst miss {100*np.max(np.abs((Z*mH + N*mn - M)/M)):.2f}% of the mass")
pr(f"   with the binding layer:    rms miss {100*np.sqrt(np.mean(rel**2)):.4f}% of the mass, worst {100*np.max(np.abs(rel)):.3f}%")
for s in ("C-12", "O-16", "Fe-56", "Pb-208", "U-238"):
    i = lab.index(s); pr(f"   {s:7s} measured {M[i]:11.2f} MeV  master table {Mp[i]:11.2f} MeV  ({Mp[i]-M[i]:+7.2f} MeV)")

head("3c. The whole chain for one atom (iron-56): every layer from the Higgs field up")
i = lab.index("Fe-56"); u_, d_ = P["up"], P["down"]
higgs_part = 26*(2*u_ + d_) + 30*(u_ + 2*d_) + 26*me
pr(f"  Higgs layer (quark and electron masses): {higgs_part:8.1f} MeV ({100*higgs_part/M[i]:.2f}% of the atom)")
pr(f"  Equilibrium layer (strong containment making 56 nucleons): {26*938.272 + 30*939.565 - (higgs_part - 26*me):8.1f} MeV")
pr(f"  Nuclear interaction layer (binding, released when iron formed): -{Bx[i]:.1f} MeV measured, -{Bp[i]:.1f} MeV from the balance formula")
pr(f"  Total: {M[i]:.1f} MeV measured. The equilibrium layer carries about {100*(1 - higgs_part/M[i]):.0f}% of the mass.")
