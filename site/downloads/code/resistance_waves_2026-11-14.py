"""Idea tested: the wave pattern is energy spreading over regions of different resistance. Three checks.
1) Resistance as refractive index n = 1 + 2GM/(r c^2) around the Sun -> light bending & Shapiro delay (Eddington; Gordon 1923 optical metric).
2) Periodic resistance -> allowed/forbidden bands (Bloch 1928; Kronig-Penney 1931), via 1D transfer matrices.
3) Random resistance -> waves trapped in place (Anderson 1958 localization)."""
import numpy as np, math
G=6.674e-11; c=2.998e8; M=1.989e30; R=6.957e8
# 1) ray through n(r): deflection = integral of grad_perp n along straight path (weak field)
b=R; z=np.linspace(-2e12,2e12,2_000_001); r=np.sqrt(b*b+z*z)
dndb=-2*G*M/(c*c)*b/r**3
defl=-np.trapezoid(dndb,z); print("1) light bending at Sun's edge from resistance n(r): %.3f arcsec (measured 1.75; Einstein 1.75; Newton-only 0.87)"%(defl*206265))
au=1.496e11; rE=au; rM=1.52*au
shap=(2*G*M/c**3)*math.log(4*rE*rM/b**2); print("   extra delay radar Earth-Mars grazing Sun (one way, Shapiro): %.0f microseconds (measured ~ 100-250 us round trip range)"%(shap*1e6))
# 2) periodic resistance: alternating layers of index n1, n2; transfer matrix trace gives bands |Tr|/2<=1 allowed
def tr(k,n1=1.0,n2=1.5,a=0.5,bb=0.5):
    def layer(n,d):
        q=k*n; return np.array([[math.cos(q*d),math.sin(q*d)/n],[-n*math.sin(q*d),math.cos(q*d)]])
    T=layer(n1,a)@layer(n2,bb); return abs(np.trace(T))/2
ks=np.linspace(0.01,12,4000); allowed=np.array([tr(k)<=1 for k in ks])
edges=np.where(np.diff(allowed.astype(int))!=0)[0]
gaps=[(round(ks[edges[i]],2),round(ks[edges[i+1]],2)) for i in range(0,len(edges)-1,2) if not allowed[edges[i]+1]]
print("2) periodic resistance -> forbidden frequency gaps (bands between them):",gaps[:5])
# 3) random resistance: transmission through N random layers decays exponentially (localization)
rng=np.random.default_rng(0)
def T_random(N,k=3.0,spread=0.3):
    T=np.eye(2)
    for _ in range(N):
        n=1+spread*rng.random(); d=0.5; q=k*n
        T=np.array([[math.cos(q*d),math.sin(q*d)/n],[-n*math.sin(q*d),math.cos(q*d)]])@T
    t=2/abs(T[0,0]+T[1,1]+1j*(T[0,1]*1-T[1,0]/1)); return t*t
for N in [10,100,1000,3000]:
    print("3) %4d random layers: transmitted fraction %.2e"%(N,np.mean([T_random(N) for _ in range(20)])))
