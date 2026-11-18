import math, numpy as np
from scipy.interpolate import PchipInterpolator
F_='font-family="Inter, Helvetica, Arial, sans-serif"'
W=2060; o=[]; BG="#f4f7fa"
def txt(x,y,s,sz=15,c="#1f2a33",w=400,a="start",halo=True,hc=None):
    h=f' paint-order="stroke" stroke="{hc or BG}" stroke-width="4"' if halo else ''
    o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{sz}" fill="{c}" font-weight="{w}" text-anchor="{a}"{h}>{s}</text>')
NU,LE,DN,UP,GO,GR,HI,RED,SF="#7a5aa6","#2f6db5","#2e8b57","#d1495b","#715800","#56616b","#444b9a","#c0392b","#1b6f8a"
nu=[3.6e-4,(3.6e-4**2+7.42e-5)**0.5,(3.6e-4**2+2.51e-3)**0.5]
le=[0.51099895e6,105.658e6,1776.93e6]; dn=[4.70e6,93.5e6,4183e6]; up=[2.16e6,1273e6,172570e6]
N=[le[i]**2/nu[i] for i in range(3)]
GRAIN,MMG,DD,H0,MAJ,DE=8.63e27,2.41e17,4.62e22,1.44e-33,8.8e3,2.24e-3
lg=math.log10
gm=lambda m: 10**(sum(lg(v) for v in m)/3)
J=lg(gm(le)/gm(nu))   # decades per turn of the colourless chain
print("jump per turn: x%.3g = %.3f decades"%(10**J,J))
print("whole chart: surface %.3g eV to grain %.3g eV = %.2f decades = %.2f turns"%(H0,GRAIN,lg(GRAIN/H0),lg(GRAIN/H0)/J))
print("middle of the whole chart sqrt(grain x surface) = %.3g eV ; neutrino middle = %.3g eV ; ratio %.2f"%((GRAIN*H0)**.5,gm(nu),gm(nu)/(GRAIN*H0)**.5))
print("turns from neutrino middle: up to grain %.3f, down to surface %.3f"%(lg(GRAIN/gm(nu))/J,lg(gm(nu)/H0)/J))
for m1 in (0.1,0.25,0.36,1.0):
    n=[m1*1e-3,((m1*1e-3)**2+7.42e-5)**.5,((m1*1e-3)**2+2.51e-3)**.5]; j=lg(gm(le)/gm(n)); print("  m1 %.2f meV: jump %.2f dec; turns up %.3f, down %.3f"%(m1,j,lg(GRAIN/gm(n))/j,lg(gm(n)/H0)/j))
steps=[("nu1>nu2",nu[0],nu[1]),("nu2>nu3",nu[1],nu[2]),("nu3>e",nu[2],le[0]),("e>mu",le[0],le[1]),("mu>tau",le[1],le[2]),("tau>N1",le[2],N[0]),("N1>N2",N[0],N[1]),("N2>N3",N[1],N[2]),("N3>grain",N[2],GRAIN)]
print("decades per third of a turn along the one spiral:",[(a,round(lg(c/b),2)) for a,b,c in steps]," (even spiral would be %.2f each)"%(J/3))
# even-spiral test: best single pitch through the 9 colourless rungs at their own directions
L=np.array([lg(v) for v in nu+le+N]); th=np.array([210+120*i for i in range(9)])/360.0
res=L-(L.mean()+J*(th-th.mean())); print("even spiral with that pitch: misses in decades:",np.round(res,2),"rms %.2f"%res.std())
# ---------- drawing ----------
CX=960; ZT=1560; YS=0.25; Y0=150+750*YS+ZT
LS,LG_=lg(H0),lg(GRAIN)
def zr(l):
    u=(l-LS)/(LG_-LS); return u*ZT, 70+u*680
def P(l,deg): z,r=zr(l); a=math.radians(deg); return CX+r*math.cos(a), Y0-z-r*math.sin(a)*YS
def ell(l,col,sw,dash="",fill="none"): z,r=zr(l); o.append(f'<ellipse cx="{CX}" cy="{Y0-z:.1f}" rx="{r:.1f}" ry="{r*YS:.1f}" fill="{fill}" stroke="{col}" stroke-width="{sw}" {dash}/>')
txt(40,50,"Each ladder's own spiral, extended without adjustment",30,w=700,halo=False)
txt(40,80,"Height is energy at true scale. Each ladder keeps its own rule: turn a third of a turn, and make the next step the same fraction of the last one. Nothing is fitted to anything else.",16,"#3b4650",halo=False)
txt(40,104,"Thick = between the three real particles. Thin = the same spiral carried on. Open dots = where further rungs would fall. Ringed dots = landings worth looking at.",16,"#3b4650",halo=False)
zt,rt=zr(LG_); zb,rb=zr(LS)
o.append(f'<polygon points="{CX-rb},{Y0} {CX-rt},{Y0-zt} {CX+rt},{Y0-zt} {CX+rb},{Y0}" fill="#eef0f4" stroke="none"/>')
ell(LG_,GO,3.5,fill="#f1eee4"); 
lm=lg(gm(nu)); zm,rm=zr(lm)
o.append(f'<polygon points="{CX-rb},{Y0} {CX-rm},{Y0-zm} {CX+rm},{Y0-zm} {CX+rb},{Y0}" fill="#e7eff3" stroke="none"/>')
o.append(f'<line x1="{CX-rb}" y1="{Y0}" x2="{CX-rt}" y2="{Y0-zt}" stroke="#9aa7b2" stroke-width="1.3"/><line x1="{CX+rb}" y1="{Y0}" x2="{CX+rt}" y2="{Y0-zt}" stroke="#9aa7b2" stroke-width="1.3"/>')
ell(LS,SF,3.5,fill="#d6e8ee")
for d in (210,330,90):
    x1,y1=P(LS,d); x2,y2=P(LG_,d); o.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#9fb3c4" stroke-width="1.1" stroke-dasharray="3 5"/>')
def rl(l,lab,col,sub="",dash='stroke-dasharray="5 4"',sw=1.8):
    ell(l,col,sw,dash+' stroke-opacity="0.75"'); x,y=P(l,200); txt(x+16,y+5,lab,14,col,700,"start",True,"#eef0f4")
    if sub: txt(x+16,y+22,sub,12.5,"#3b4650",400,"start",True,"#eef0f4")
rl(LG_,"THE GRAIN  8.6 × 10¹⁸ GeV",GO,"the top","",3.5)
rl(lg(DD),"D (coloured far-side parts) 4.6 × 10¹³ GeV",GR)
rl(lg(MMG),"MMG 2.4 × 10⁸ GeV",RED)
rl(lg(1e11),"W, Z, Higgs",GO)
rl(lg(MAJ),"majoron 8.8 keV",NU)
rl(lm,"halfway line (3.5 meV); dark energy's scale (2.2 meV) is here too","#1f2a33","",'stroke-dasharray="2 4"',1.6)
ell(LS,SF,3.5); txt(CX-rb-18,Y0+4,"THE SURFACE  1.4 × 10⁻³³ eV",14.5,SF,700,"end"); txt(CX-rb-18,Y0+22,"the bottom",12.5,"#3b4650",400,"end")
# each ladder's own spiral, extended without adjustment: L(t) = Linf - C k^t, one rung every third of a turn
def law(m):
    L=[lg(v) for v in m]; d1,d2=L[1]-L[0],L[2]-L[1]; k=d2/d1; Li=L[2]+d2*k/(1-k); return k,Li,Li-L[0]
NOTE={}
def strand(t0,t1,fn,col,sw,extra=""):
    pts=[P(fn(t),210+120*t) for t in np.arange(t0,t1+1e-9,0.01)]
    o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="{col}" stroke-width="{sw}" stroke-linejoin="round" stroke-linecap="round" {extra}/>')
for name,m,col in (("neutrinos",nu,NU),("leptons",le,LE),("hidden",N,HI),("up",up,UP),("down",dn,DN)):
    k,Li,C=law(m); fn=lambda t,Li=Li,C=C,k=k: min(max(Li-C*k**t,LS),LG_)
    if k<1: t0=math.log((Li-LS)/C)/math.log(k); t1=2+math.log(0.02/(C*k**2*(1-k)))/math.log(k)
    else: t1=math.log((Li-LG_)/C)/math.log(k); t0=math.log(0.02/(abs(C)*(k-1)))/math.log(k)
    strand(t0,0,fn,col,2.0,'stroke-opacity="0.7"'); strand(2,t1,fn,col,2.0,'stroke-opacity="0.7"'); strand(0,2,fn,col,4.6)
    for t in range(int(math.floor(t0))-1,int(math.ceil(t1))+2):
        if t in (0,1,2): continue
        l=Li-C*k**t
        if LS-0.5<=l<=LG_+0.5 and abs(C*k**t*(1-k))>0.5:
            x,y=P(min(max(l,LS),LG_),210+120*t); o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5.5" fill="#fff" stroke="{col}" stroke-width="2.2"/>'); NOTE[(name,t)]=(x,y,l)
    print(name,"k %.3f"%k,"ceiling" if k<1 else "floor","%.3g eV"%10**Li)
    z,r=zr(Li); x=CX+r+26; y=Y0-z
    o.append(f'<line x1="{x-14:.1f}" y1="{y:.1f}" x2="{x+6:.1f}" y2="{y:.1f}" stroke="{col}" stroke-width="3"/>')
    lab={"neutrinos":"neutrino spiral: ceiling 0.45 eV","leptons":"lepton spiral: ceiling 42 GeV","hidden":"hidden-partner spiral: ceiling 4.2 × 10¹² GeV","up":"up-quark spiral: ceiling 2.3 × 10⁹ GeV","down":"down-quark spiral: FLOOR 76 eV (it runs the other way)"}[name]
    txt(x+12,y+5,lab,14,col,700)
def call(key,text,dx,dy,a="start"):
    if key in NOTE:
        x,y,l=NOTE[key]; o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="11" fill="none" stroke="#1f2a33" stroke-width="2.4"/>'); txt(x+dx,y+dy,text,14,"#1f2a33",700,a)
call(("neutrinos",-4),"neutrino spiral, 4 steps down: lands on the surface",16,-14)
call(("hidden",-2),"hidden spiral, 2 steps down: neutrino 1's energy (0.33 meV), a third of a turn round",16,22)
call(("down",7),"down spiral, 7th step: 3.6 × 10¹⁸ GeV, just under the grain",-16,26,"end")
call(("down",6),"down spiral, 6th step: 1.2 × 10¹³ GeV, near D",16,-12)
call(("up",-2),"up spiral, 2 steps down: 11 meV, beside neutrino 2",16,20)
def dot(v,th,col,name,dx,dy,a,r=7.5,hollow=False):
    x,y=P(lg(v),th); o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{"#fff" if hollow else col}" stroke="{col if hollow else "#fff"}" stroke-width="2"/>'); txt(x+dx,y+dy,name,13.5,col,700,a)
A3=[210,330,450]
LBL=[(-14,6,"end"),(14,6,"start"),(14,-8,"start")]
for i in range(3):
    dx,dy,a=LBL[i]
    dot(nu[i],A3[i],NU,["neutrino 1  0.36 meV (the theory's value)","neutrino 2  8.6 meV","neutrino 3  50 meV"][i],dx,dy,a,6.5)
    q=[(-14,20,"end"),(14,24,"start"),(14,16,"start")][i]; dot(le[i],A3[i],LE,["electron 0.511 MeV","muon 105.7 MeV","tau 1.777 GeV"][i],*q)
    q=[(-14,-12,"end"),(14,8,"start"),(-14,-2,"end")][i]; dot(dn[i],A3[i],DN,["down 4.7 MeV","strange 93.5 MeV","bottom 4.18 GeV"][i],*q)
    q=[(-14,4,"end"),(14,-6,"start"),(14,-8,"start")][i]; dot(up[i],A3[i],UP,["up 2.2 MeV","charm 1.27 GeV","top 172.6 GeV"][i],*q)
    dot(N[i],A3[i],HI,["hidden 1  7 × 10⁵ GeV","hidden 2  1.3 × 10⁹ GeV","hidden 3  6.3 × 10¹⁰ GeV"][i],dx,dy,a,8.5)
ty=Y0+rb*YS+50; tb=len(o); o.append(""); BG="#ffffff"
txt(40,ty+36,"What the extended spirals show",21,w=700,halo=False)
Ls=["Four spirals (neutrinos, leptons, hidden partners, up quarks) tighten as they climb. Carried on, each winds up to a ceiling and never reaches the grain. Carried down, each opens out and reaches the surface in a few steps.",
"The down-quark spiral runs the other way: it has a floor (76 eV) and it is the only one that climbs all the way to the grain. Its 6th and 7th steps fall near D and just under the grain.",
"Two sharp landings, with the extra part-turn you expected: the neutrino spiral reaches the surface in exactly 4 steps if the lightest neutrino is 0.363 meV, and the hidden spiral comes down onto",
"neutrino 1's energy in 2 steps if it is 0.358 meV. The theory's own value, from the neutrino rule, is 0.36 meV. Other landings (near D, the grain, neutrino 2) are no better than chance.",
"The open dots are not particles: a fourth charged lepton at 8 GeV, for example, is ruled out. So the spirals connect as waves, at most; they do not add rungs."]
for i,l in enumerate(Ls): txt(40,ty+70+i*27,l,15.5,halo=False)
H=int(ty+70+5*27+20)
o[tb]=f'<rect x="20" y="{ty}" width="{W-40}" height="{H-ty-20}" fill="#ffffff" rx="10" stroke="#d5dde4"/>'
out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" {F_}><rect width="{W}" height="{H}" fill="#f4f7fa"/>']+o+['</svg>']
open("multi.svg","w").write("\n".join(out)); print("H",H)
