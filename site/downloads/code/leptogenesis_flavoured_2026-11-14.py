"""Flavoured leptogenesis for the 3x3-link hidden partners (14 Nov 2026).
Regimes: 1e9 < T < 1e12 GeV -> two flavours (tau, and a = e+mu); T < 1e9 -> three flavours.
  N3 (6.3e10) and N2 (1.3e9) act in the two-flavour regime; N1 (7e5) in the three-flavour regime.
Flavoured CP asymmetries (Covi, Roulet & Vissani 1996; Nardi, Nir, Roulet & Racker 2006; Abada et al. 2006):
  eps_ia = -(1/(8 pi (hh+)_ii)) sum_{j!=i} { Im[h_ia h*_ja (hh+)_ij] (f_V+f_S)(x) + Im[h_ia h*_ja (hh+)_ji] /(1-x) },  x = M_j^2/M_i^2
Washout: K_ia = K_i |h_ia|^2/(hh+)_ii; efficiency kappa(K) (Buchmuller, Di Bari & Plumacher 2005).
N3's asymmetry in flavour a is split into the part N2 sees (projection p) and the part it can't see (Barbieri et al. 2000;
Engelhard, Grossman, Nardi & Nir 2007; Di Bari 2005; Bertuzzo, Di Bari & Marzola 2011). Flavour-coupling matrix ignored (O(1)).
eta_B ~ 0.96e-2 x total lepton asymmetry. Observed 6.1e-10 (Planck 2018). Mixing: NuFIT 6.0; delta and Majorana phases scanned."""
import math, numpy as np
rng = np.random.default_rng(29)
v = 174.0; mstar = 1.08e-12
m = np.array([0.00036, 0.0087, 0.0501])*1e-9
ml = np.array([0.51099895e-3, 0.1056584, 1.77693]); M = ml**2/m
def fV(x):
    if x > 1e4: return -math.sqrt(x)*(1/(2*x) - 1/(6*x*x))
    return math.sqrt(x)*(1 - (1+x)*math.log1p(1/x))
def fS(x): return math.sqrt(x)/(1-x)
def kappa(K):
    if K < 1e-12: return 0.0
    zB = 2 + 4*K**0.13*math.exp(-2.5/K); return 2/(zB*K)*(1 - math.exp(-K*zB/2))
def rot(i, j, z):
    R = np.eye(3, dtype=complex); c, s = np.cos(z), np.sin(z); R[i, i] = c; R[j, j] = c; R[i, j] = s; R[j, i] = -s; return R
def pmns(s12, s13, s23, d, a21, a31):
    c12, c13, c23 = math.sqrt(1-s12), math.sqrt(1-s13), math.sqrt(1-s23); s12, s13, s23 = map(math.sqrt, (s12, s13, s23))
    e = np.exp(-1j*d)
    U = np.array([[c12*c13, s12*c13, s13*e],
                  [-s12*c23 - c12*s23*s13/e, c12*c23 - s12*s23*s13/e, s23*c13],
                  [s12*s23 - c12*c23*s13/e, -c12*s23 - s12*c23*s13/e, c23*c13]])
    return U@np.diag([1, np.exp(1j*a21/2), np.exp(1j*a31/2)])
def eta(R, U):
    h = (np.sqrt(M)[:, None]*R*np.sqrt(m)[None, :])@U.conj().T/v      # h[i, alpha]
    hh = h@h.conj().T
    K = np.real(np.diag(hh))*v**2/(M*mstar)
    eps = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if j == i: continue
            x = M[j]**2/M[i]**2
            for a in range(3):
                t1 = np.imag(h[i, a]*np.conj(h[j, a])*hh[i, j])*(fV(x)+fS(x))
                t2 = np.imag(h[i, a]*np.conj(h[j, a])*hh[j, i])/(1-x)
                eps[i, a] += -(t1+t2)/(8*math.pi*np.real(hh[i, i]))
    P = np.abs(h)**2/np.real(np.diag(hh))[:, None]; Kf = K[:, None]*P
    # N3 in two flavours: tau (2), a = e+mu (0,1)
    e3t, e3a = eps[2, 2], eps[2, 0]+eps[2, 1]
    K3t, K3a = Kf[2, 2], Kf[2, 0]+Kf[2, 1]
    Nt = e3t*kappa(K3t); Na = e3a*kappa(K3a)
    # N2 washout: tau directly; a split into seen (p) and unseen (1-p)
    K2t, K2a = Kf[1, 2], Kf[1, 0]+Kf[1, 1]
    l3 = h[2, :2]/np.linalg.norm(h[2, :2]); l2 = h[1, :2]/np.linalg.norm(h[1, :2])
    p = abs(np.vdot(l2, l3))**2
    Nt = Nt*math.exp(-3*math.pi*K2t/8)
    Na_seen = Na*p*math.exp(-3*math.pi*K2a/8); Na_unseen = Na*(1-p)
    # N2 own production
    Nt += eps[1, 2]*kappa(K2t); Na_seen += (eps[1, 0]+eps[1, 1])*kappa(K2a)
    # below 1e9: three flavours. Split the a-asymmetry into e and mu by N2's direction (seen part) and the orthogonal (unseen part)
    w2 = np.abs(l2)**2; wperp = 1 - w2
    Ne = Na_seen*w2[0] + Na_unseen*wperp[0]; Nm = Na_seen*w2[1] + Na_unseen*wperp[1]
    # N1 washout in three flavours
    Ne *= math.exp(-3*math.pi*Kf[0, 0]/8); Nm *= math.exp(-3*math.pi*Kf[0, 1]/8); Nt *= math.exp(-3*math.pi*Kf[0, 2]/8)
    return 0.96e-2*(Ne+Nm+Nt), K, p
results = []
for s23 in (0.47, 0.561):
    for y in (0.3, 1.0):
        vals = []
        for _ in range(30000):
            zs = rng.uniform(0, 2*math.pi, 3) + 1j*rng.uniform(-y, y, 3)
            R = rot(1, 2, zs[0])@rot(0, 2, zs[1])@rot(0, 1, zs[2])
            U = pmns(0.307, 0.02195, s23, rng.uniform(0, 2*math.pi), rng.uniform(0, 2*math.pi), rng.uniform(0, 2*math.pi))
            e_, K, p = eta(R, U); vals.append(abs(e_))
        vals = np.array(vals)
        print(f"sin^2 th23 = {s23}, misalignment |Im| <= {y}: max eta_B = {vals.max():.2e}, median {np.median(vals):.1e}, "
              f"fraction >= 6.1e-10: {np.mean(vals >= 6.1e-10):.4f}")
        results.append(vals.max())
print(f"Best flavoured result {max(results):.2e} vs observed 6.1e-10 (ratio {6.1e-10/max(results):.0f}x short)" if max(results) < 6.1e-10
      else f"Flavoured result reaches the observed value: best {max(results):.2e}")
print("Unflavoured best: 4.6e-12.")

# Local search from the best random points (is the true maximum higher than random sampling found?)
def params_eta(pv, s23):
    zs = pv[0:3] + 1j*pv[3:6]
    R = rot(1, 2, zs[0])@rot(0, 2, zs[1])@rot(0, 1, zs[2])
    U = pmns(0.307, 0.02195, s23, pv[6], pv[7], pv[8])
    return abs(eta(R, U)[0])
for s23 in (0.47, 0.561):
    starts = []
    for _ in range(4000):
        pv = np.concatenate([rng.uniform(0, 2*math.pi, 3), rng.uniform(-0.5, 0.5, 3), rng.uniform(0, 2*math.pi, 3)])
        starts.append((params_eta(pv, s23), pv))
    starts.sort(key=lambda t: -t[0])
    best_val = 0; best_pv = None
    for val, pv in starts[:8]:
        step = 0.3
        for it in range(1500):
            trial = pv + rng.normal(scale=step, size=9); trial[3:6] = np.clip(trial[3:6], -1, 1)
            tv = params_eta(trial, s23)
            if tv > val: val, pv = tv, trial
            if it % 300 == 299: step *= 0.5
        if val > best_val: best_val, best_pv = val, pv
    print(f"local search, sin^2 th23 = {s23}: max eta_B = {best_val:.2e} (observed 6.1e-10); "
          f"misalignment imaginary parts {np.round(best_pv[3:6],2)}, delta = {math.degrees(best_pv[6])%360:.0f} deg")
