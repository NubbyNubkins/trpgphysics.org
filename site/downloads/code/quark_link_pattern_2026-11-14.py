"""Quark sectors in the same relational (circulant) form as the leptons (14 Nov 2026).
sqrt-mass eigenvalues M0 (1 + 2 r cos(delta + 2 pi k/3)) give Koide Q = (1 + 2 r^2)/3, so r^2 = (3Q - 1)/2.
Leptons: Q = 2/3 -> r^2 = 1/2. Quarks: Q at a common scale (QCD running leaves Q unchanged; quark_common_scale script).
Uncertainties: PDG 2024 mass errors, Monte Carlo. Look-elsewhere: how often would simple sixths match this well by chance?"""
import numpy as np, math
rng = np.random.default_rng(0); N = 200000
# PDG 2024 MS-bar values (GeV) and 1-sigma errors, run to 2 GeV-equivalent ratios are scale-free for Q, so use one consistent set:
# values at mu = MZ from the quark_common_scale script (same running for all, so relative errors carry over)
mZ = {"u": (0.00128, 0.00128*0.07/2.16*1.0*10), "d": (0.00278, 0.00278*0.07/4.70*10),
      "s": (0.0552, 0.0552*0.008/0.0935*10/1.0), "c": (0.658, 0.658*0.0046/1.273), "b": (2.889, 2.889*0.007/4.183), "t": (170.3, 170.3*0.7/162.5)}
# relative errors (PDG 2024): u +0.49/-0.26 MeV of 2.16 (~17%), d +0.48/-0.17 of 4.70 (~7%), s +0.8/-0.3 of 93.5 (~6%),
# c 4.6/1273 (0.4%), b 7/4183 (0.2%), t(MS-bar) ~0.7/162.5 (0.4%)
rel = {"u": 0.17, "d": 0.07, "s": 0.06, "c": 0.0036, "b": 0.0017, "t": 0.0043}
cen = {k: v[0] for k, v in mZ.items()}
def Q(a, b, c): return (a+b+c)/(np.sqrt(a)+np.sqrt(b)+np.sqrt(c))**2
def r2(q): return (3*q-1)/2
def delta(ms):
    s = np.sqrt(ms); M = s.mean(); return math.acos((s[-1]/M-1)/(2*math.sqrt(r2(Q(*ms)))))
samp = {k: cen[k]*(1+rel[k]*rng.standard_normal(N)) for k in cen}
for k in samp: samp[k] = np.abs(samp[k])
qu = Q(samp["u"], samp["c"], samp["t"]); qd = Q(samp["d"], samp["s"], samp["b"])
for name, q, target, label in [("up-type (u,c,t)", qu, 5/6, "5/6"), ("down-type (d,s,b)", qd, 2/3, "2/3")]:
    rr = r2(q)
    print(f"{name}: Q = {np.median(q):.4f} +/- {q.std():.4f};  r^2 = {np.median(rr):.4f} +/- {rr.std():.4f};  "
          f"simple value {label} = {target:.4f} is {(np.median(rr)-target)/rr.std():+.1f} sigma away")
print("leptons: r^2 = 1/2 (Koide). Tested idea: r^2 = 3/6, 4/6, 5/6 for leptons, down-type, up-type -> FAILS (both quark sectors are ~7 sigma off).")
print("twist angles: leptons 2/9 = %.4f; up-type %.4f; down-type %.4f" % (2/9, delta(np.array([cen['u'], cen['c'], cen['t']])), delta(np.array([cen['d'], cen['s'], cen['b']]))))
# Look-elsewhere: chance that two random values (uniform on the allowed range r^2 in [0.5, 1)) both land within
# their measured distance of SOME sixth or ninth or twelfth (simple fractions a researcher might try)
fr = sorted({a/b for b in (2, 3, 4, 6, 8, 9, 12) for a in range(1, b) if 0.5 <= a/b < 1})
du = min(abs(np.median(r2(qu))-f) for f in fr); dd = min(abs(np.median(r2(qd))-f) for f in fr)
x = rng.uniform(0.5, 1.0, (N, 2)); near = lambda v, d: np.min(np.abs(v[:, None]-np.array(fr)[None, :]), axis=1) <= d
print(f"look-elsewhere: {len(fr)} simple fractions tried; chance both random values land as close: {np.mean(near(x[:,0],du)&near(x[:,1],dd)):.3f}")
