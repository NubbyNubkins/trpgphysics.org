"""Can the hidden right-handed partners (an earlier test, 3x3 link) make the matter excess? (14 Nov 2026)
Unflavoured estimate with standard approximations:
  Casas-Ibarra (2001): h = (1/v) sqrt(M) R sqrt(m) U^dagger, R complex orthogonal ("misalignment" between the link and the partners)
  CP asymmetry (Covi, Roulet & Vissani 1996): eps_i = -(1/8pi) sum_j Im[(hh+)_ij^2]/(hh+)_ii [f_V + f_S](M_j^2/M_i^2)
  decay parameter K_i = (hh+)_ii v^2 / (M_i m*), m* = 1.08e-3 eV
  efficiency kappa(K) (Buchmuller, Di Bari & Plumacher 2005); washout by lighter partners exp(-3 pi K_j / 8) (Di Bari 2005)
  eta_B ~ 0.96e-2 * eps * kappa ; observed eta_B = 6.1e-10 (Planck 2018)
Credit: Sakharov; Fukugita & Yanagida (1986); Casas & Ibarra; Covi, Roulet & Vissani; Buchmuller, Di Bari & Plumacher; Di Bari;
Davidson & Ibarra; Nardi, Nir, Roulet & Racker (2006) and Abada et al. (2006) (flavour effects)."""
import math, numpy as np
rng = np.random.default_rng(28)
v = 174.0; mstar = 1.08e-12      # GeV
m = np.array([0.00036, 0.0087, 0.0501])*1e-9
ml = np.array([0.51099895e-3, 0.1056584, 1.77693])
M = ml**2/m
def fV(x):
    if x > 1e4: return -math.sqrt(x)*(1/(2*x) - 1/(6*x*x))      # series: avoids round-off
    return math.sqrt(x)*(1 - (1+x)*math.log1p(1/x))
def fS(x): return math.sqrt(x)/(1-x)
def kappa(K):
    zB = 2 + 4*K**0.13*math.exp(-2.5/K)
    return 2/(zB*K)*(1 - math.exp(-K*zB/2))
def rot(i, j, z):
    R = np.eye(3, dtype=complex); c, s = np.cos(z), np.sin(z)
    R[i, i] = c; R[j, j] = c; R[i, j] = s; R[j, i] = -s; return R
def run(R):
    hh = (np.sqrt(M)[:, None]*(R*m)@R.conj().T*np.sqrt(M)[None, :])/v**2
    K = np.real(np.diag(hh))*v**2/(M*mstar)
    eps = np.zeros(3)
    for i in range(3):
        for j in range(3):
            if j == i: continue
            x = M[j]**2/M[i]**2
            eps[i] += -np.imag(hh[i, j]**2)/np.real(hh[i, i])*(fV(x)+fS(x))/(8*math.pi)
    eta = np.zeros(3)
    for i in range(3):
        wash = math.prod(math.exp(-3*math.pi*K[j]/8) for j in range(i))   # lighter partners wash out
        eta[i] = 0.96e-2*eps[i]*kappa(K[i])*wash
    return eta, eps, K
eta0, eps0, K0 = run(np.eye(3))
print(f"Partners M = {M[0]:.2e}, {M[1]:.2e}, {M[2]:.2e} GeV")
print(f"Perfectly aligned (R = 1, exactly the earlier assumption): CP asymmetries = {eps0} -> NO matter excess")
print(f"  decay parameters K = {np.round(K0,2)} (K>1 = strong washout)")
print(f"  washout of partner 3's excess by partners 1,2: {math.exp(-3*math.pi*(K0[0]+K0[1])/8):.1e}")
best = (0, None)
for y in (0.3, 1.0):
    vals = []
    for _ in range(40000):
        zs = rng.uniform(0, 2*math.pi, 3) + 1j*rng.uniform(-y, y, 3)
        R = rot(1, 2, zs[0])@rot(0, 2, zs[1])@rot(0, 1, zs[2])
        eta, eps, K = run(R)
        tot = abs(eta.sum()); vals.append(tot)
        if tot > best[0]: best = (tot, (eta, eps, K, zs))
    vals = np.array(vals)
    print(f"Misalignment allowed up to imaginary part {y}: max |eta_B| = {vals.max():.1e}; fraction reaching 6.1e-10: {np.mean(vals>6.1e-10):.4f}")
eta, eps, K, zs = best[1]
print(f"  best case: eta from partners 1,2,3 = {eta}, K = {np.round(K,1)}")
print("Observed eta_B = 6.1e-10.")
# cross-check against the Davidson-Ibarra closed form for partner 1
R = rot(1,2,zs[0])@rot(0,2,zs[1])@rot(0,1,zs[2])
DI = -3*M[0]/(16*math.pi*v**2)*np.imag(np.sum(m**2*R[0]**2))/np.sum(m*abs(R[0])**2)
print(f"  check: partner-1 asymmetry {eps[0]:.3e} vs Davidson-Ibarra closed form {DI:.3e}")
