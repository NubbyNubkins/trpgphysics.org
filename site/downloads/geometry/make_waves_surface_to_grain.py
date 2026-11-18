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
txt(40,50,"The waves stretched at the surface and squished at the grain, with antimatter",30,w=700,halo=False)
txt(40,80,"Height is energy. Width is the length of one wave: widest at the surface, where one wave is as long as the whole surface, narrowest at the grain, where one wave is one grain long.",16,"#3b4650",halo=False)
txt(40,104,"Solid = matter. Dashed with open dots = antimatter: the same wave at the same height, wound the other way. Grey = six whole turns fitted exactly between the two ends.",16,"#3b4650",halo=False)
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
LEV=["the surface","nothing known here","nothing known here","3.5 meV: the neutrinos (middle 5.4 meV) and dark energy (2.2 meV)","47.5 MeV: middle of the leptons is 45.8 MeV, 4% away","6.4 × 10⁸ GeV: hidden partners (middle 3.9 × 10⁸) and the MMG (2.4 × 10⁸)","the grain"]
for k in range(7):
    l=LS+k*J; z,r=zr(l)
    if 0<k<6: ell(l,"#1f2a33",1.3,'stroke-dasharray="2 5"',op=0.7)
    x=CX+r+22; y=Y0-z; o.append(f'<line x1="{x-12:.1f}" y1="{y:.1f}" x2="{x+8:.1f}" y2="{y:.1f}" stroke="#1f2a33" stroke-width="1.8"/>')
    txt(x+16,y+1,("turn %d: "%k if 0<k<6 else ("THE GRAIN  8.6 × 10¹⁸ GeV" if k==6 else "THE SURFACE  1.4 × 10⁻³³ eV"))+(LEV[k] if 0<k<6 else ""),15,[SF,"#6b7884","#6b7884",NU,LE,HI,GO][k],700)
    txt(x+16,y+20,"one wave here is "+wl(10**l)+" long",13,"#3b4650")
ell(LG_,GO,3.5,fill="#f1eee4")
# six whole turns, matter and antimatter
def helix(sign,col,sw,dash=""):
    pts=[]
    for i in range(0,2161,3):
        l=LS+(LG_-LS)*i/2160; pts.append(P(l,90+sign*i))
    o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="{col}" stroke-width="{sw}" stroke-opacity="0.55" {dash}/>')
helix(+1,"#6b7884",1.8); helix(-1,"#6b7884",1.4,'stroke-dasharray="5 5"')
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
txt(CX-RW+10,Y0+RW*YS+30,"Open dots: the antimatter partners (positron, anti-muon, anti-quarks). On the back line the two sit on the same spot. Neutrinos and hidden partners are their own partners.",14,"#3b4650",halo=False)
# ---- unrolled strip ----
sy=Y0+RW*YS+150; x0,x1=120,W-120
txt(40,sy-70,"The same thing unrolled: one wave from the surface to the grain",22,w=700,halo=False)
txt(40,sy-46,"Six turns. Each turn is really 13.5 billion times shorter than the one before; it is drawn far gentler so you can see it. Blue = matter, red dashed = antimatter (the mirror wave).",15,"#3b4650",halo=False)
q=0.74; tot=sum(q**k for k in range(6)); xs=[x0]
for k in range(6): xs.append(xs[-1]+(x1-x0)*q**k/tot)
for sign,col,extra in ((1,LE,""),(-1,RED,'stroke-dasharray="8 6"')):
    pts=[]
    for k in range(6):
        for i in range(0,101):
            x=xs[k]+(xs[k+1]-xs[k])*i/100; amp=70*q**(k+i/100); pts.append((x,sy+70-sign*amp*math.sin(2*math.pi*i/100)))
    o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="{col}" stroke-width="3" {extra}/>')
o.append(f'<line x1="{x0}" y1="{sy+70}" x2="{x1}" y2="{sy+70}" stroke="#1f2a33" stroke-width="1.2"/>')
SL=["surface","turn 1","turn 2","turn 3: neutrinos","turn 4: leptons, quarks","turn 5: hidden partners","grain"]
for k in range(7):
    o.append(f'<line x1="{xs[k]:.1f}" y1="{sy-10}" x2="{xs[k]:.1f}" y2="{sy+150}" stroke="#9aa7b2" stroke-width="1.2" stroke-dasharray="3 4"/>')
    txt(xs[k],sy+170,SL[k],13.5,[SF,"#6b7884","#6b7884",NU,LE,HI,GO][k],700,"middle" if 0<k<6 else ("start" if k==0 else "end"),halo=False)
txt(x0,sy+192,"stretched",14,SF,600,halo=False); txt(x1,sy+192,"squished",14,GO,600,"end",halo=False)
ty=sy+220; tb=len(o); o.append(""); BG="#ffffff"
txt(40,ty+36,"What comes out",21,w=700,halo=False)
Ls=["Stretched and squished is literal. A wave's length is set by its energy: at the surface one wave is as long as the whole surface; at the grain one wave is one grain long. Everything else sits in between.",
"A value does come out. If exactly six turns fit between the two ends, the fourth level is 47.5 MeV. The middle of the three leptons is 45.8 MeV: 4% away, with nothing put in but the surface and the grain.",
"The same six turns put level 3 at 3.5 meV, where the neutrinos and dark energy are, and level 5 at 6.4 × 10⁸ GeV, where the hidden partners and the MMG are (each within a factor of about 2).",
"Caution: I tried 3 to 8 turns and six is the one that fits, so the chance of luck is about 2%, not 0.3%. It is the same relation Weinberg noticed in 1972 between the cosmos and particle sizes. A lead, not a pass.",
"Tension: six exact turns want the lightest neutrino near 0.10 meV; the last picture's landings want 0.36 meV. Both cannot be exact.",
"Antimatter is the same wave at the same height, turning the other way. When the two meet, the turning (the charge) cancels to zero and the energy is left over as light. That is measured: electron + positron gives two flashes of 0.511 MeV."]
for i,l in enumerate(Ls): txt(40,ty+70+i*27,l,15.5,halo=False)
H=int(ty+70+6*27+20)
o[tb]=f'<rect x="20" y="{ty}" width="{W-40}" height="{H-ty-20}" fill="#ffffff" rx="10" stroke="#d5dde4"/>'
out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" {F_}><rect width="{W}" height="{H}" fill="#f4f7fa"/>']+o+['</svg>']
open("funnel.svg","w").write("\n".join(out)); print("H",H)
