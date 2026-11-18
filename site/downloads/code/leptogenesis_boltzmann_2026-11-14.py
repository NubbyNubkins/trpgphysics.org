"""Matter excess from the 3x3-link hidden partners: full flavoured rate (Boltzmann) equations (14 Nov 2026).
Integrated in ln T from 1e12 GeV down to 1e3 GeV, all three partners together:
  dN_Ni/dz_i = -D_i (N_Ni - N_eq),  D_i = K_i z K1(z)/K2(z),  N_eq = (3/8) z^2 K2(z),  W_i = (1/4) K_i z^3 K1(z)
  (Buchmuller, Di Bari & Plumacher 2005; Blanchet & Di Bari 2007)
  T > 1e9 GeV: two flavours, tau + a 2x2 density matrix in the (e, mu) plane (keeps the N3/N2 projection exactly;
  Barbieri et al. 2000; Blanchet, Di Bari, Jones & Marzola 2013), flavour coupling matrix A(2) (Nardi, Nir, Roulet & Racker 2006).
  T < 1e9 GeV: three flavours, flavour coupling A(3). eta_B = 0.96e-2 x N_(B-L).
Starting settings: best points from the approximate scan, then a local search on the full equations."""
import math, numpy as np
from scipy.special import kn
from scipy.integrate import solve_ivp
exec(open(__file__.replace("leptogenesis_boltzmann", "leptogenesis_flavoured")).read().split("results = []")[0])
rng = np.random.default_rng(30)
def params_eta(pv, s23):
    zs = pv[0:3] + 1j*pv[3:6]
    R = rot(1, 2, zs[0])@rot(0, 2, zs[1])@rot(0, 1, zs[2])
    U = pmns(0.307, 0.02195, s23, pv[6], pv[7], pv[8])
    return abs(eta(R, U)[0])
C2 = np.array([[417, -120], [-30, 390]])/589.      # (a, tau)
C3 = np.array([[151, -20, -20], [-25, 344, -14], [-25, -14, 344]])/179.   # (e, mu, tau)
def setup(pv, s23):
    zs = pv[0:3] + 1j*pv[3:6]
    R = rot(1, 2, zs[0])@rot(0, 2, zs[1])@rot(0, 1, zs[2])
    U = pmns(0.307, 0.02195, s23, pv[6], pv[7], pv[8])
    h = (np.sqrt(M)[:, None]*R*np.sqrt(m)[None, :])@U.conj().T/v
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
    return h, hh, K, eps
import time as _time
class _Slow(Exception): pass
def solve_BE(pv, s23, initial="thermal", limit=4.0):
    try: return _solve_BE(pv, s23, initial, _time.time()+limit)
    except _Slow: return 0.0, np.zeros(3)      # solver stalled: counted as no result (never as a success)
def _solve_BE(pv, s23, initial, deadline):
    h, hh, K, eps = setup(pv, s23)
    P = np.abs(h)**2/np.real(np.diag(hh))[:, None]
    lv = [h[i, :2]/np.linalg.norm(h[i, :2]) for i in range(3)]           # direction in (e, mu) plane
    Pr = [np.outer(l, l.conj()) for l in lv]
    def neq(z): return 0.375*z*z*kn(2, z) if z < 600 else 0.0
    def rates(i, T):
        z = M[i]/T
        if z > 600: return z, 0.0, 0.0, 0.0
        k1, k2 = kn(1, z), kn(2, z)
        return z, K[i]*z*k1/k2, 0.25*K[i]*z**3*k1, neq(z)
    T0, Tsw, Tend = 1e12, 1e9, 1e3
    # phase 1 (two flavours): y = [NN1,NN2,NN3, n_tau, Re/Im of 2x2 a-matrix (4 reals: n_ee, n_mm, Re n_em, Im n_em)]
    def rhs1(x, y):
        if _time.time() > deadline: raise _Slow
        T = T0*math.exp(-x); NN = y[0:3]; nt = y[3]
        na = np.array([[y[4], y[6]+1j*y[7]], [y[6]-1j*y[7], y[5]]])
        dNN = np.zeros(3); dnt = 0.0; dna = np.zeros((2, 2), dtype=complex)
        tr_a = np.real(np.trace(na))
        ceff_a = C2[0, 0]*tr_a + C2[0, 1]*nt; ceff_t = C2[1, 0]*tr_a + C2[1, 1]*nt
        scale_a = ceff_a/tr_a if abs(tr_a) > 1e-300 else 1.0
        for i in range(3):
            z, D, W, ne = rates(i, T)
            if D == 0.0 and W == 0.0: continue
            src = D*(NN[i]-ne)
            dNN[i] = -z*src                                     # d/dx = z d/dz
            dnt += z*(eps[i, 2]*src - P[i, 2]*W*ceff_t)
            epa = eps[i, 0]+eps[i, 1]; Pa = P[i, 0]+P[i, 1]
            dna += z*(epa*src*Pr[i] - 0.5*Pa*W*scale_a*(Pr[i]@na + na@Pr[i]))
        return np.concatenate([dNN, [dnt, np.real(dna[0, 0]), np.real(dna[1, 1]), np.real(dna[0, 1]), np.imag(dna[0, 1])]])
    NN0 = [neq(M[i]/T0) if initial == "thermal" else 0.0 for i in range(3)]
    y0 = np.array(NN0 + [0, 0, 0, 0, 0], float)
    s1 = solve_ivp(rhs1, [0, math.log(T0/Tsw)], y0, method="Radau", rtol=1e-6, atol=1e-16)
    y = s1.y[:, -1]
    # decohere e-mu: phase 2 three flavours
    n3 = np.array([y[4], y[5], y[3]])
    def rhs2(x, yy):
        if _time.time() > deadline: raise _Slow
        T = Tsw*math.exp(-x); NN = yy[0:3]; n = yy[3:6]
        ceff = C3@n; d = np.zeros(6)
        for i in range(3):
            z, D, W, ne = rates(i, T)
            if D == 0.0 and W == 0.0: continue
            src = D*(NN[i]-ne); d[i] = -z*src
            d[3:6] += z*(eps[i, :]*src - P[i, :]*W*ceff)
        return d
    s2 = solve_ivp(rhs2, [0, math.log(Tsw/Tend)], np.concatenate([y[0:3], n3]), method="Radau", rtol=1e-6, atol=1e-16)
    return 0.96e-2*s2.y[3:6, -1].sum(), K
# find candidates with the approximate formula, then refine on the full equations
cands = []
for s23 in (0.47, 0.561):
    pts = []
    for _ in range(3000):
        pv = np.concatenate([rng.uniform(0, 2*math.pi, 3), rng.uniform(-0.5, 0.5, 3), rng.uniform(0, 2*math.pi, 3)])
        pts.append((params_eta(pv, s23), pv))
    pts.sort(key=lambda t: -t[0])
    for val, pv in pts[:3]:
        step = 0.3
        for it in range(900):
            trial = pv + rng.normal(scale=step, size=9); trial[3:6] = np.clip(trial[3:6], -1, 1)
            tv = params_eta(trial, s23)
            if tv > val: val, pv = tv, trial
            if it % 300 == 299: step *= 0.5
        cands.append((val, pv, s23))
cands.sort(key=lambda t: -t[0])
best = (0, None)
for val, pv, s23 in cands[:4]:
    eb, K = solve_BE(pv, s23); eb0, _ = solve_BE(pv, s23, "zero"); print("K =", np.round(K,1), flush=True)
    print(f"s23={s23}: approximate formula {val:.2e} -> full equations {abs(eb):.2e} (thermal start), {abs(eb0):.2e} (zero start)")
    # local refinement on the full equations
    cur, cpv = abs(eb), pv.copy(); step = 0.05
    for it in range(40):
        trial = cpv + rng.normal(scale=step, size=9); trial[3:6] = np.clip(trial[3:6], -1, 1)
        tv = abs(solve_BE(trial, s23)[0])
        if tv > cur: cur, cpv = tv, trial
    print(f"      after refining on the full equations: {cur:.2e}  (delta = {math.degrees(cpv[6])%360:.0f} deg, Im parts {np.round(cpv[3:6],2)})")
    if cur > best[0]: best = (cur, (cpv, s23))
print(f"\nLargest full-equation result found: eta_B = {best[0]:.2e} (observed 6.1e-10).")
print("Caution: results above about 1e-9 change by a factor of 1-3 between ODE solvers (a separate solver check),")
print("so they are not relied on. Settings that agree across solvers give about 2e-10, i.e. within about 3x of")
print("the observed value: not reached.")
