# Extend each ladder's own spiral without adjustment (turn a third, step by the ladder's own shrink factor k)
import numpy as np
lg=np.log10
def nus(m1): return [m1*1e-3,((m1*1e-3)**2+7.42e-5)**.5,((m1*1e-3)**2+2.51e-3)**.5]
le=[0.51099895e6,105.6584e6,1776.93e6]; up=[2.16e6,1273e6,172570e6]; dn=[4.70e6,93.5e6,4183e6]
GRAIN,MMG,DD,H0,MAJ=8.63e27,2.41e17,4.62e22,1.44e-33,8.8e3
LS,LG=lg(H0),lg(GRAIN)
def law(m):
    L=lg(np.array(m)); d1,d2=L[1]-L[0],L[2]-L[1]; k=d2/d1; Linf=L[2]+d2*k/(1-k); C=Linf-L[0]; return L,k,Linf,C   # L(t)=Linf-C k^t, t=0,1,2 the rungs
def run(m1=0.36,show=True):
    nu=nus(m1); N=[le[i]**2/nu[i] for i in range(3)]
    lad={"neutrinos":nu,"leptons":le,"hidden":N,"up":up,"down":dn}
    marks={"surface":(LS,None),"nu1":(lg(nu[0]),210),"nu2":(lg(nu[1]),330),"nu3":(lg(nu[2]),90),"majoron":(lg(MAJ),None),"e":(lg(le[0]),210),"u":(lg(up[0]),210),"d":(lg(dn[0]),210),"mu":(lg(le[1]),330),"s":(lg(dn[1]),330),"c":(lg(up[1]),330),"tau":(lg(le[2]),90),"b":(lg(dn[2]),90),"W/Z/H":(lg(1e11),None),"t":(lg(up[2]),90),
           "hidden1":(lg(N[0]),210),"MMG":(lg(MMG),None),"hidden2":(lg(N[1]),330),"hidden3":(lg(N[2]),90),"D":(lg(DD),None),"grain":(LG,None)}
    out={}
    for n,m in lad.items():
        L,k,Linf,C=law(m); pts=[]
        for t in list(range(-1,-12,-1))+list(range(3,14)):
            l=Linf-C*k**t
            if LS-0.6<=l<=LG+0.6 and abs(C*k**t*(1-k))>0.15: pts.append((t,l,(210+120*t)%360))
        out[n]=(L,k,Linf,C,sorted(pts))
        if show:
            print(f"\n{n}: shrink k={k:.3f}; {'ceiling' if k<1 else 'floor'} at 10^{Linf:.2f} eV = {10**Linf:.3g} eV")
            for t,l,a in sorted(pts):
                near=sorted(((abs(l-v[0]),nm,v[1]) for nm,v in marks.items()))[0]
                print(f"   step {t:+d}: {10**l:.3g} eV at direction {a:3.0f}   nearest landmark {near[1]} (x{10**near[0]:.2f}{', same direction' if near[2]==a else (', any direction' if near[2] is None else ', other direction')})")
    return out,marks,lad
out,marks,lad=run()
def hits(out,marks,own,tol=0.4):
    n=0
    for nm,(L,k,Linf,C,pts) in out.items():
        for t,l,a in pts:
            for mk,(v,d) in marks.items():
                if mk in own.get(nm,()): continue
                if abs(l-v)<tol and (d is None or d==a): n+=1; break
    return n
own={"neutrinos":("nu1","nu2","nu3"),"leptons":("e","mu","tau"),"hidden":("hidden1","hidden2","hidden3"),"up":("u","c","t"),"down":("d","s","b")}
h=hits(out,marks,own); tot=sum(len(v[4]) for v in out.values())
print(f"\nLandings within x2.5 of a landmark (right direction where the landmark has one): {h} of {tot} extension points")
# chance: move each ladder's middle rung at random (keeps ends), recount
rng=np.random.default_rng(3); cnt=[]
for _ in range(4000):
    o2={}
    for nm,m in lad.items():
        L=lg(np.array(m)); f=rng.uniform(0.25,0.75); Lm=L[0]+f*(L[2]-L[0]); d1,d2=Lm-L[0],L[2]-Lm; k=d2/d1
        if abs(k-1)<1e-3: k=1.001
        Linf=L[2]+d2*k/(1-k); C=Linf-L[0]; pts=[]
        for t in list(range(-1,-12,-1))+list(range(3,14)):
            l=Linf-C*k**t
            if LS-0.6<=l<=LG+0.6 and abs(C*k**t*(1-k))>0.15: pts.append((t,l,(210+120*t)%360))
        o2[nm]=(L,k,Linf,C,pts)
    cnt.append((hits(o2,marks,own),sum(len(v[4]) for v in o2.values())))
cnt=np.array(cnt); print(f"random middle rungs: mean {cnt[:,0].mean():.2f} landings of {cnt[:,1].mean():.1f} points; chance of {h} or more = {np.mean(cnt[:,0]>=h)*100:.0f}%")
print("\nDependence on the lightest neutrino: neutrino step -4 vs surface; hidden step -2 vs neutrino 1")
for m1 in (0.1,0.25,0.36,0.5,1.0):
    o3,mk,_=run(m1,False); Ln=o3["neutrinos"]; Lh=o3["hidden"]
    print(f"  m1 {m1:4.2f}: neutrino spiral 4 steps down = 10^{Ln[2]-Ln[3]*Ln[1]**-4:.2f} eV (surface 10^{LS:.2f}); hidden 2 steps down = {10**(Lh[2]-Lh[3]*Lh[1]**-2)*1e3:.3g} meV (neutrino 1 = {m1})")
