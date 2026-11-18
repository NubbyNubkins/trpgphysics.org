"""Run all six quark masses to one common energy (2-loop QCD), then retest the Koide balance angle. 14 Nov 2026.
Inputs PDG 2024: alpha_s(MZ)=0.1180; u,d,s MS-bar at 2 GeV; c(mc)=1.273; b(mb)=4.183; t MS-bar(mt) ~ 162.5 GeV (from pole 172.57).
Running: d a/dln mu^2 = -b0 a^2 - b1 a^3 ; d ln m/d ln mu^2 = -g0 a - g1 a^2 ; a = alpha_s/pi (Tarrach; Chetyrkin; Vermaseren)."""
import math
from scipy.integrate import solve_ivp
MZ=91.1876
def b(nf): return (11-2*nf/3)/4, (102-38*nf/3)/16
def g(nf): return 1.0, (202/3-20*nf/9)/16
def nf_at(mu): return 6 if mu>162.5 else 5 if mu>4.183 else 4 if mu>1.273 else 3
def run(a0,mu0,mu1,lnm=0.0):
    # piecewise across thresholds, continuous matching (adequate at 2 loops for this test)
    th=[1.273,4.183,162.5]; pts=sorted([mu0,mu1]+[t for t in th if min(mu0,mu1)<t<max(mu0,mu1)],reverse=(mu1<mu0))
    a,l=a0,lnm
    for x,y in zip(pts,pts[1:]):
        nf=nf_at(math.sqrt(x*y)); b0,b1=b(nf); g0,g1=g(nf)
        f=lambda t,s:[-b0*s[0]**2-b1*s[0]**3, -g0*s[0]-g1*s[0]**2]
        s=solve_ivp(f,[math.log(x*x),math.log(y*y)],[a,l],rtol=1e-10,atol=1e-12).y[:,-1]; a,l=s
    return a,l
aZ=0.1180/math.pi
def a_at(mu): return run(aZ,MZ,mu)[0]
def m_to(m0,mu0,mu1): return m0*math.exp(run(a_at(mu0),mu0,mu1)[1])
inp={"u":(2.16e-3,2.0),"d":(4.70e-3,2.0),"s":(0.0935,2.0),"c":(1.273,1.273),"b":(4.183,4.183),"t":(162.5,162.5)}
print("alpha_s(2 GeV)=%.4f  alpha_s(mt)=%.4f"%(a_at(2.0)*math.pi,a_at(162.5)*math.pi))
def Qang(v):
    Q=sum(v)/sum(math.sqrt(x) for x in v)**2; return Q, math.degrees(math.acos(math.sqrt(1/(3*Q))))
for mu in (2.0,MZ,162.5,1000.0):
    m={k:m_to(m0,mu0,mu) for k,(m0,mu0) in inp.items()}
    up=Qang([m['u'],m['c'],m['t']]); dn=Qang([m['d'],m['s'],m['b']])
    print("at %7.1f GeV: u %.5f d %.5f s %.4f c %.3f b %.3f t %.1f GeV | up-type Q=%.4f angle %.2f | down-type Q=%.4f angle %.2f"%(mu,m['u'],m['d'],m['s'],m['c'],m['b'],m['t'],up[0],up[1],dn[0],dn[1]))
