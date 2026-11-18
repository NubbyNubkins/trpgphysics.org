"""The 'reverse engineer by scaling backwards'. The elementary masses (Higgs links) and the three
strengths are run from the top's energy up to the Planck energy with the Standard Model's running equations.
Test 1: does any quark ladder become a 45-degree cone (like the leptons) at some energy, and is that energy one the theory already has?
Test 2: run the top's link DOWN from many starting values at high energy: does it get pulled to its measured value (the fixed-point pull)?
Credit: Gross & Wilczek 1973, Politzer 1973 (running); Machacek & Vaughn 1983-85; Arason et al. 1992 (Yukawa running); Pendleton & Ross 1981,
Hill 1981 (fixed point); Buttazzo et al. 2013 (inputs at the top's energy); Xing, Zhang & Zhou 2008/2012 (masses at the top's energy);
Koide 1981; Xing & Zhang 2006 (the lepton 45 degrees in running masses)."""
import math, numpy as np
from scipy.integrate import solve_ivp
def pr(*a): print(*a, flush=True)
def head(t): pr("\n" + "="*104); pr(t); pr("="*104)
mt0 = 172.5; v = 246.22; k = 16*math.pi**2
g0 = [0.4626, 0.6478, 1.1666]                                     # g1 (GUT normalised), g2, g3 at the top's energy (Buttazzo et al. 2013)
mu_ = np.array([1.242, 627.8, 0.0]); md_ = np.array([2.678, 53.36, 2873.0]); me_ = np.array([0.48657, 102.718, 1746.24])   # MeV
yu = np.sqrt(2)*mu_/1000/v; yu[2] = 0.9369; yd = np.sqrt(2)*md_/1000/v; ye = np.sqrt(2)*me_/1000/v
b1 = np.array([41/10, -19/6, -7.0]); B2 = np.array([[199/50, 27/10, 44/5], [9/10, 35/6, 12], [11/10, 9/2, -26]]); ct = np.array([17/10, 3/2, 2])
def rhs(t, y):
    g = y[:3]; u = y[3:6]; d = y[6:9]; e = y[9:12]; g2 = g**2
    T = 3*np.sum(u**2) + 3*np.sum(d**2) + np.sum(e**2)
    dg = g**3*b1/k + g**3/k**2*(B2 @ g2 - ct*u[2]**2)
    du = u/k*(1.5*(u**2 - d**2) + T - (17/20*g2[0] + 9/4*g2[1] + 8*g2[2]))
    dd = d/k*(1.5*(d**2 - u**2) + T - (1/4*g2[0] + 9/4*g2[1] + 8*g2[2]))
    de = e/k*(1.5*e**2 + T - (9/4*g2[0] + 9/4*g2[1]))
    return [*dg, *du, *dd, *de]
tmax = math.log(1.22e19/mt0)
sol = solve_ivp(rhs, [0, tmax], [*g0, *yu, *yd, *ye], rtol=1e-10, atol=1e-14, dense_output=True)
def Q(m): r = np.sqrt(np.abs(m)); return np.sum(np.abs(m))/np.sum(r)**2
def ang(q): return math.degrees(math.acos(math.sqrt(1/(3*q))))
head("1. The three ladders' cone openings from the top's energy to the Planck energy (45 degrees = the lepton rule, Q = 2/3)")
pr(f"  {'energy (GeV)':>14s} {'leptons':>9s} {'down':>9s} {'up':>9s}   strengths 1/a1 1/a2 1/a3")
marks = {"top's energy": mt0, "far-side pairs": 2.1e8, "hidden partners": 6.3e10, "strong = weak": None, "Planck": 1.2e19}
E = np.logspace(math.log10(mt0), 19, 400); rows = []
for En in E:
    y = sol.sol(math.log(En/mt0)); u, d, e = y[3:6], y[6:9], y[9:12]
    rows.append((En, ang(Q(e)), ang(Q(d)), ang(Q(u)), *(4*math.pi/y[:3]**2)))
rows = np.array(rows)
i23 = np.argmin(np.abs(rows[:, 5] - rows[:, 6])); marks["strong = weak"] = rows[i23, 0]
for En in (mt0, 1e3, 1e5, 2.1e8, 6.3e10, rows[i23, 0], 1e16, 1.2e19):
    y = sol.sol(math.log(En/mt0)); u, d, e = y[3:6], y[6:9], y[9:12]
    pr(f"  {En:14.3e} {ang(Q(e)):9.3f} {ang(Q(d)):9.3f} {ang(Q(u)):9.3f}   {4*math.pi/y[0]**2:6.1f} {4*math.pi/y[1]**2:6.1f} {4*math.pi/y[2]**2:6.1f}")
pr(f"  Range over the whole climb: leptons {rows[:,1].min():.3f} to {rows[:,1].max():.3f}; down {rows[:,2].min():.3f} to {rows[:,2].max():.3f}; up {rows[:,3].min():.3f} to {rows[:,3].max():.3f} degrees.")
for name, col in (("down", 2), ("up", 3)):
    hit = np.where(np.diff(np.sign(rows[:, col] - 45)))[0]
    pr(f"  {name} ladder reaches 45 degrees: " + (", ".join(f"{rows[h,0]:.2e} GeV" for h in hit) if len(hit) else "never, between the top's energy and the Planck energy"))
pr("  Why: the strong force pulls all quark masses down by almost the same factor as the energy rises, so ratios inside a ladder barely move.")
pr("  The opening of a ladder depends only on its mass ratios, so it hardly changes with energy. Scaling up does not turn quarks into cones.")
y = sol.sol(math.log(2.1e8/mt0)); pr(f"  At the far-side scale: b/tau = {y[8]/y[11]:.3f}, s/mu = {y[7]/y[10]:.3f}, d/e = {y[6]/y[9]:.3f}; top link {y[5]:.3f}")
y = sol.sol(tmax); pr(f"  At the Planck energy:  b/tau = {y[8]/y[11]:.3f}, s/mu = {y[7]/y[10]:.3f}, d/e = {y[6]/y[9]:.3f}; top link {y[5]:.3f}")

head("2. The top's link run DOWN from many starting values (the fixed-point pull)")
def run_down(yt_hi, E_hi):
    t_hi = math.log(E_hi/mt0); ys = sol.sol(t_hi).copy(); ys[5] = yt_hi
    s = solve_ivp(rhs, [t_hi, 0], ys, rtol=1e-10, atol=1e-14)
    return s.y[5, -1]
for E_hi, name in ((rows[i23, 0], "strong = weak"), (1e16, "1e16 GeV"), (1.2e19, "Planck")):
    out = [(y0, run_down(y0, E_hi)) for y0 in (0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0)]
    pr(f"  start at {name:14s}: " + ", ".join(f"{a:g} -> {b:.3f}" for a, b in out))
    pr(f"      a large start lands at about {out[-1][1]:.3f} (top mass {out[-1][1]*v/math.sqrt(2):.0f} GeV in this scheme); measured 0.937 ({0.937*v/math.sqrt(2):.0f} GeV)")
pr("  Reading: any large top link at high energy is pulled down toward one value, so the top's link near 1 needs no fine choice. But the pull")
pr("  lands 36-41% above the measured top's link in the Standard Model alone (222-229 GeV against 163 GeV in this scheme; the known result of")
pr("  Hill 1981). The measured top needs a start of only about 0.41-0.5 at high energy, which is not on the pull. A 'full glass' link of 1")
pr("  at high energy lands at about 1.19-1.20 (about 207 GeV): also too heavy.")
