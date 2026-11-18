import math, numpy as np
F_='font-family="Inter, Helvetica, Arial, sans-serif"'
W=2060; o=[]; BG="#f4f7fa"
def txt(x,y,s,sz=15,c="#1f2a33",w=400,a="start",halo=True,hc=None):
    h=f' paint-order="stroke" stroke="{hc or BG}" stroke-width="4"' if halo else ''
    o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{sz}" fill="{c}" font-weight="{w}" text-anchor="{a}"{h}>{s}</text>')
NU,LE,DN,UP,GO,GR,HI,RED,SF="#7a5aa6","#2f6db5","#2e8b57","#d1495b","#715800","#56616b","#444b9a","#c0392b","#1b6f8a"
lg=math.log10
nu=[3.6e-4,(3.6e-4**2+7.42e-5)**0.5,(3.6e-4**2+2.51e-3)**0.5]
le=[0.51099895e6,105.658e6,1776.93e6]; dn=[4.70e6,93.5e6,4183e6]; up=[2.16e6,1273e6,172570e6]
N=[le[i]**2/nu[i] for i in range(3)]
GRAIN,H0=8.63e27,1.44e-33; LS,LG_=lg(H0),lg(GRAIN); J=(LG_-LS)/6
CX=820; ZT=1500; YS=0.24; RW=640; RT=46; Y0=150+ZT
def zr(l): u=(l-LS)/(LG_-LS); return u*ZT, RT+(1-u)*(RW-RT)
def P(l,deg): z,r=zr(l); a=math.radians(deg); return CX+r*math.cos(a), Y0-z-r*math.sin(a)*YS
def ell(l,col,sw,dash="",fill="none",op=1): z,r=zr(l); o.append(f'<ellipse cx="{CX}" cy="{Y0-z:.1f}" rx="{r:.1f}" ry="{r*YS:.1f}" fill="{fill}" stroke="{col}" stroke-width="{sw}" stroke-opacity="{op}" {dash}/>')
txt(40,50,"Each spiral carried on from its three particles, across the whole funnel",30,w=700,halo=False)
txt(40,80,"Height is energy. Width is the length of one wave: widest at the surface, where one wave is as long as the whole surface, narrowest at the grain, where one wave is one grain long.",16,"#3b4650",halo=False)
txt(40,104,"Nothing is adjusted: each spiral keeps turning a third of a turn per step, and each step is the same fraction of the one before, as set by its own three particles.",16,"#3b4650",halo=False)
o.append(f'<polygon points="{CX-RW},{Y0} {CX-RT},{Y0-ZT} {CX+RT},{Y0-ZT} {CX+RW},{Y0}" fill="#eef0f4" stroke="none"/>')
ell(LS,SF,3.5,fill="#d6e8ee")
o.append(f'<line x1="{CX-RW}" y1="{Y0}" x2="{CX-RT}" y2="{Y0-ZT}" stroke="#9aa7b2" stroke-width="1.3"/><line x1="{CX+RW}" y1="{Y0}" x2="{CX+RT}" y2="{Y0-ZT}" stroke="#9aa7b2" stroke-width="1.3"/>')
for d in (210,330,90):
    x1,y1=P(LS,d); x2,y2=P(LG_,d); o.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#9fb3c4" stroke-width="1" stroke-dasharray="3 5"/>')
def wl(E):
    m=1.2398e-6/E
    for u_,f_ in (("billion light-years",9.461e24),("light-years",9.461e15),("km",1e3),("m",1),("mm",1e-3),("thousandths of a mm",1e-6),("nm",1e-9),("fm",1e-15)):
        if m>=f_:
            v=m/f_; return (("{:,.0f}".format(float("%.2g"%v)) if v>=100 else "%.2g"%v))+" "+u_
    e=int(math.floor(math.log10(m))); return ("%.1f × 10"%(m/10**e))+str(e).translate(str.maketrans("-0123456789","⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))+" m"
for k in range(7):
    l=LS+k*J; z,r=zr(l)
    if 0<k<6: ell(l,"#1f2a33",1.1,'stroke-dasharray="2 6"',op=0.5)
    x=CX-r-22; y=Y0-z; txt(x,y+5,["THE SURFACE  1.4 × 10⁻³³ eV","level 1","level 2","level 3: 3.5 meV","level 4: 47.5 MeV","level 5: 6.4 × 10⁸ GeV","THE GRAIN  8.6 × 10¹⁸ GeV"][k],14,[SF,"#6b7884","#6b7884","#6b7884","#6b7884","#6b7884",GO][k],700 if k in (0,6) else 500,"end") if k else txt(CX+r+22,y+5,"THE SURFACE  1.4 × 10⁻³³ eV",14,SF,700)
ell(LG_,GO,3.5,fill="#f1eee4")
# six whole turns, matter and antimatter
def helix(sign,col,sw,dash=""):
    pts=[]
    for i in range(0,2161,3):
        l=LS+(LG_-LS)*i/2160; pts.append(P(l,90+sign*i))
    o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="{col}" stroke-width="{sw}" stroke-opacity="0.55" {dash}/>')

# each ladder's own spiral carried on from its three points: L(t) = Linf - C k^t
def law(m):
    L=[lg(v) for v in m]; d1,d2=L[1]-L[0],L[2]-L[1]; k=d2/d1; Li=L[2]+d2*k/(1-k); return k,Li,Li-L[0]
CEIL=[]
for name,m,col in (("hidden partners",N,HI),("up quarks",up,UP),("down quarks",dn,DN),("leptons",le,LE),("neutrinos",nu,NU)):
    k,Li,C=law(m); fn=lambda t,Li=Li,C=C,k=k: min(max(Li-C*k**t,LS),LG_)
    if k<1: t0=math.log((Li-LS)/C)/math.log(k); t1=2+math.log(0.03/(C*k**2*(1-k)))/math.log(k)
    else: t1=math.log((Li-LG_)/C)/math.log(k); t0=math.log(0.03/(abs(C)*(k-1)))/math.log(k)
    for sign,extra,sw in ((+1,'stroke-opacity="0.85"',2.2),(-1,'stroke-dasharray="6 5" stroke-opacity="0.55"',1.5)):
        pts=[P(fn(t),(210+120*t) if sign>0 else (180-(210+120*t))) for t in np.arange(t0,t1+1e-9,0.01)]
        o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="{col}" stroke-width="{sw}" stroke-linejoin="round" {extra}/>')
    ell(Li,col,2.6,op=0.9); CEIL.append((name,Li,col,k))
    # where it meets the end
    te=t0 if k<1 else t1; x,y=P(LS if k<1 else LG_,210+120*te); o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="8" fill="{col}" stroke="#1f2a33" stroke-width="2"/>')
# ladders: matter thick, antimatter mirrored thin
def ladder(m,col,names,sw=4.2,lab=None):
    L=[lg(v) for v in m]; d1,d2=L[1]-L[0],L[2]-L[1]
    f=lambda t: L[0]+d1*t+(d2-d1)*t*(t-1)/2
    for sign,extra,w_ in ((+1,"",sw),(-1,'stroke-dasharray="7 5"',sw-2.0)):
        pts=[P(f(t),(210+120*t) if sign>0 else (180-(210+120*t))) for t in np.arange(0,2.001,0.02)]
        o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="{col}" stroke-width="{w_}" stroke-linecap="round" {extra}/>')
    for i in range(3):
        x,y=P(L[i],210+120*i); xa,ya=P(L[i],180-(210+120*i))
        if i<2: o.append(f'<circle cx="{xa:.1f}" cy="{ya:.1f}" r="6" fill="#fff" stroke="{col}" stroke-width="2.4"/>')
        else: o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="11" fill="none" stroke="{col}" stroke-width="2.2" stroke-dasharray="3 3"/>')
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{col}" stroke="#fff" stroke-width="1.8"/>')
        dx,dy,a=(lab or [(-14,6,"end"),(14,6,"start"),(14,-10,"start")])[i]
        if names[i]: txt(x+dx,y+dy,names[i],13.5,col,700,a)
ladder(N,HI,["hidden 1","hidden 2","hidden 3  6.3 × 10¹⁰ GeV"])
ladder(up,UP,["up","charm","top 172.6 GeV"],3.4,[(-14,-6,"end"),(14,-8,"start"),(14,-10,"start")])
ladder(dn,DN,["down","strange","bottom"],3.4,[(-14,-20,"end"),(14,8,"start"),(-14,-4,"end")])
ladder(le,LE,["electron 0.511 MeV","muon 105.7 MeV","tau 1.777 GeV"],4.2,[(-14,20,"end"),(14,24,"start"),(14,16,"start")])
ladder(nu,NU,["neutrino 1","neutrino 2","neutrino 3  50 meV"])
CL={"neutrinos":"neutrino spiral: comes off the surface, stops at its ceiling, 0.45 eV","leptons":"lepton spiral: comes off the surface, stops at 42 GeV","up quarks":"up-quark spiral: comes off the surface, stops at 2.3 × 10⁹ GeV","hidden partners":"hidden-partner spiral: comes off the surface, stops at 4.2 × 10¹² GeV","down quarks":"down-quark spiral: comes down from the grain, stops at its floor, 76 eV"}
used=[]
for name,Li,col,k in sorted(CEIL,key=lambda c:-c[1]):
    z,r=zr(Li); x=CX+r+(150 if name=='hidden partners' else 22); y=Y0-z
    while any(abs(y-u)<24 for u in used): y+=24
    used.append(y)
    o.append(f'<line x1="{CX+r:.1f}" y1="{Y0-z:.1f}" x2="{x+10:.1f}" y2="{y:.1f}" stroke="{col}" stroke-width="2"/>'); txt(x+16,y+5,CL[name],15,col,700)
txt(CX-RW+10,Y0+RW*YS+30,"Thick = between the three real particles. Thin = the same spiral carried on by its own rule. Dashed = antimatter, the same spiral wound the other way. Big dots = where a spiral meets the surface or the grain.",14,"#3b4650",halo=False)
ty=Y0+RW*YS+60; tb=len(o); o.append(""); BG="#ffffff"
txt(40,ty+36,"What it looks like",21,w=700,halo=False)
Ls=["No single spiral runs all the way from top to bottom. Carried on by its own rule, each one reaches one end and stops short of the other.",
"Four come off the surface: neutrinos, leptons, up quarks and hidden partners. Each climbs through its three particles, tightens, and winds up into a ring, its ceiling. None reaches the grain.",
"One comes down from the grain: the down quarks. It tightens the other way and winds down into a ring at 76 eV. It never reaches the surface.",
"So the picture is waves coming off the surface and a wave coming to it from the grain: both. They pass through each other above 76 eV. All the charged particles sit there; the neutrinos sit below it.",
"Near the surface the spirals are wide and loose (stretched). At each ceiling they pile up into a tight ring (squished). The squish happens at each spiral's own ceiling, not at the grain.",
"Antimatter is each spiral's mirror twin, wound the other way, crossing its partner twice every turn."]
for i,l in enumerate(Ls): txt(40,ty+70+i*27,l,15.5,halo=False)
H=int(ty+70+6*27+20)
o[tb]=f'<rect x="20" y="{ty}" width="{W-40}" height="{H-ty-20}" fill="#ffffff" rx="10" stroke="#d5dde4"/>'
out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" {F_}><rect width="{W}" height="{H}" fill="#f4f7fa"/>']+o+['</svg>']
open("funnel2.svg","w").write("\n".join(out)); print("H",H)
