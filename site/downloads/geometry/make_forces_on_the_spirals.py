import math, numpy as np
W=2400; H=2190; o=[]; BG="#f4f7fa"
def txt(x,y,s,sz=15,c="#1f2a33",w=400,a="start",halo=True):
    h=f' paint-order="stroke" stroke="{BG}" stroke-width="4"' if halo else ''
    o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{sz}" fill="{c}" font-weight="{w}" text-anchor="{a}"{h}>{s}</text>')
NU,LE,DN,UP,GO,GR,HI,RED,SF,LT="#7a5aa6","#2f6db5","#2e8b57","#d1495b","#715800","#56616b","#444b9a","#c0392b","#1b6f8a","#c98a00"
EL,WK,ST,GV="#c98a00","#7a5aa6","#c0392b","#1f2a33"
lg=math.log10
GRAIN,H0=8.633e27,1.438e-33; LS,LG_=lg(H0),lg(GRAIN)
CX=700; ZT=1560; YS=0.22; RW=430; RT=38; Y0=260+ZT
def zr(l): u=(l-LS)/(LG_-LS); return u*ZT, RT+(1-u)*(RW-RT)
def Y(l): return Y0-zr(l)[0]
def P(l,deg): z,r=zr(l); a=math.radians(deg); return CX+r*math.cos(a), Y0-z-r*math.sin(a)*YS
def ell(l,col,sw,dash="",fill="none",op=1): z,r=zr(l); o.append(f'<ellipse cx="{CX}" cy="{Y0-z:.1f}" rx="{r:.1f}" ry="{r*YS:.1f}" fill="{fill}" stroke="{col}" stroke-width="{sw}" stroke-opacity="{op}" {dash}/>')
txt(40,52,"Where the forces are on the spirals",30,w=700,halo=False)
txt(40,84,"The spirals are the things. The forces are not on a spiral: they are what joins the spirals, what rings them, and what the cone itself is made of.",16,"#3b4650",halo=False)
txt(40,108,"Left: the cone with the five spirals and each force drawn where it acts. Right: how strong each force is at every height of the same cone.",16,"#3b4650",halo=False)
o.append(f'<polygon points="{CX-RW},{Y0} {CX-RT},{Y0-ZT} {CX+RT},{Y0-ZT} {CX+RW},{Y0}" fill="#eef0f4"/>')
ell(LS,SF,3.5,fill="#d6e8ee")
o.append(f'<line x1="{CX-RW}" y1="{Y0}" x2="{CX-RT}" y2="{Y0-ZT}" stroke="{GV}" stroke-width="2.6"/><line x1="{CX+RW}" y1="{Y0}" x2="{CX+RT}" y2="{Y0-ZT}" stroke="{GV}" stroke-width="2.6"/>')
o.append(f'<line x1="{CX}" y1="{Y0}" x2="{CX}" y2="{Y0-ZT}" stroke="#9aa7b2" stroke-width="1" stroke-dasharray="3 6"/>')
nu=[3.6e-4,8.62e-3,5.01e-2]; le=[0.51099895e6,105.658e6,1776.93e6]; dn=[4.70e6,93.5e6,4183e6]; up=[2.16e6,1273e6,172570e6]; N=[le[i]**2/nu[i] for i in range(3)]
# binding rings (push and pull forces)
for E,col,name in ((3.73e3,EL,"atom ring: electric pull"),(2.35e8,ST,"proton ring: strong pull")):
    ell(lg(E),col,4.5,op=0.9)
# far-side levels
for E,col in ((2.41e17,HI),(4.62e22,HI)): ell(lg(E),col,1.6,'stroke-dasharray="2 4"')
ell(LG_,GO,3.5,fill="#f1eee4")
def ladder(m,col,off=0):
    L=[lg(v) for v in m]; d1,d2=L[1]-L[0],L[2]-L[1]; f=lambda t: L[0]+d1*t+(d2-d1)*t*(t-1)/2
    pts=[P(f(t),210+off+120*t) for t in np.arange(0,2.001,0.02)]
    o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="{col}" stroke-width="4" stroke-linecap="round"/>')
    return [P(L[i],210+off+120*i) for i in range(3)]
pN=ladder(N,HI); pU=ladder(up,UP,8); pD=ladder(dn,DN,-8); pL=ladder(le,LE); pV=ladder(nu,NU)
# weak force: the hop between partner spirals
for a,b in list(zip(pU,pD))+list(zip(pL,pV)):
    o.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" stroke="{WK}" stroke-width="2.4" stroke-dasharray="7 5"/>')
for pts,col,names in ((pN,HI,("hidden 1","hidden 2","hidden 3")),(pU,UP,("up","charm","top")),(pD,DN,("down","strange","bottom")),(pL,LE,("electron","muon","tau")),(pV,NU,("ν1","ν2","ν3"))):
    for (x,y),n in zip(pts,names):
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6.5" fill="{col}" stroke="#fff" stroke-width="1.6"/>')
        txt(x+(11 if x>=CX else -11),y+(17 if col in (DN,LE) else -9),n,12.5,col,700,'start' if x>=CX else 'end')
# carriers on the centre line
for E,col,n in ((80.4e9,WK,"W"),(91.2e9,WK,"Z"),(125.2e9,GR,"Higgs")):
    o.append(f'<rect x="{CX-6}" y="{Y(lg(E))-6:.1f}" width="12" height="12" fill="{col}" stroke="#fff" stroke-width="1.5" transform="rotate(45 {CX} {Y(lg(E)):.1f})"/>')
o.append(f'<rect x="{CX-7}" y="{Y0-7}" width="14" height="14" fill="{EL}" stroke="#fff" stroke-width="1.5" transform="rotate(45 {CX} {Y0})"/>')
# left labels with leaders
XL=40
SLOT=[0]
def lab(l,title,sub,col,dy=0,xr=None):
    y=Y0-ZT+10+SLOT[0]*118; SLOT[0]+=1; x1=(CX-zr(l)[1]) if xr is None else xr
    o.append(f'<polyline points="{XL+250},{y-5:.1f} {x1-8:.1f},{Y(l):.1f}" fill="none" stroke="{col}" stroke-width="1.2" stroke-opacity="0.7"/>')
    txt(XL,y-8,title,15.5,col,700); 
    for i,sline in enumerate(sub): txt(XL,y+10+i*17,sline,13,"#3b4650")
lab(LG_,"GRAVITY: the cone itself",["the slope from surface to grain;","its strength is set by the grain"],GV,dy=30)
lab(lg(4.62e22),"where the three strengths meet",["coloured far-side parts (D)"],HI,dy=10)
lab(lg(2.41e17),"MMG level",["the weak line bends here"],HI,dy=40)
lab(lg(1.6e11),"W, Z and Higgs: on the centre line",["weak carriers 80 and 91 GeV;","the link that gives weight, 125 GeV;","the surface's resistance switches on"],WK,dy=-30,xr=CX)
lab(lg(2.35e8),"STRONG: the proton ring",["rings the first quark rungs; locks","at 155 MeV; never reaches past the ring"],ST,dy=10)
lab(lg(5e5),"WEAK: the hop between spirals",["dashed purple: up to down, and","electron-kind to neutrino. It joins","two spirals; it is a change of view"],WK,dy=70,xr=(pL[0][0]+pV[0][0])/2-40)
lab(lg(3.73e3),"ELECTRIC: the atom ring",["rings the electron and the nucleus;","endless reach, carried by light"],EL,dy=60)
SLOT[0]=11.2
lab(LS,"LIGHT and GLUONS: at the surface",["the carriers with no weight;","light is the surface moving"],EL,dy=-40,xr=CX)
txt(CX,Y0-ZT-RT*YS-14,"the grain",14,GO,700,"middle"); txt(CX,Y0+RW*YS+30,"the surface",15,SF,700,"middle")
# right: strengths against height
GX0,GX1=1320,2040
def X(v): return GX0+(GX1-GX0)*(v+45)/45     # v = log10(strength), -45..0
o.append(f'<rect x="{GX0}" y="{Y0-ZT}" width="{GX1-GX0}" height="{ZT}" fill="#fff" stroke="#c9d3dc"/>')
for v in range(-45,1,5):
    o.append(f'<line x1="{X(v):.1f}" y1="{Y0-ZT}" x2="{X(v):.1f}" y2="{Y0}" stroke="#e3e8ee" stroke-width="1"/>'); txt(X(v),Y0+22,"1" if v==0 else f"10^{v}".replace("^-","⁻").translate(str.maketrans("0123456789","⁰¹²³⁴⁵⁶⁷⁸⁹")).replace("1⁰","10") if False else ("1" if v==0 else "1/10"+str(-v).translate(str.maketrans("0123456789","⁰¹²³⁴⁵⁶⁷⁸⁹"))),12.5,"#56616b",a="middle",halo=False)
txt((GX0+GX1)/2,Y0+48,"strength of each force at that height (1 = full strength)",14.5,"#3b4650",600,"middle",halo=False)
txt(GX0,Y0-ZT-16,"HOW STRONG EACH FORCE IS, ALL THE WAY UP THE CONE",13,"#6b7884",700,halo=False)
MZ=91.19e9; L0=math.log
def gauge(E):
    t=L0(E/MZ); a1,a2,a3=59.0,29.6,8.45; b=(4.1,-19/6,-7.0); M=2.41e17
    if E<=M: return [a1-b[0]/(2*math.pi)*t, a2-b[1]/(2*math.pi)*t, a3-b[2]/(2*math.pi)*t]
    tm=L0(M/MZ); t2=L0(E/M); b2=(5.3,-19/6+2,-7.0)
    return [a1-b[0]/(2*math.pi)*tm-b2[0]/(2*math.pi)*t2, a2-b[1]/(2*math.pi)*tm-b2[1]/(2*math.pi)*t2, a3-b[2]/(2*math.pi)*tm-b2[2]/(2*math.pi)*t2]
def line(pts,col,sw=3,dash=""): o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="{col}" stroke-width="{sw}" stroke-linecap="round" {dash}/>')
Es=np.logspace(lg(MZ),lg(4.62e22),200)
for i,col in ((0,EL),(1,WK),(2,ST)): line([(X(-lg(gauge(E)[i])),Y(lg(E))) for E in Es],col)
# below the Z: electric nearly flat to the surface; strong grows and locks; weak hidden (short reach)
line([(X(lg(1/137.036)),Y0),(X(lg(1/137.036)),Y(lg(0.511e6))),(X(lg(1/127.9)),Y(lg(MZ)))],EL)
def as_(E):   # strong strength below the Z, one loop, 5 then 4 then 3 kinds
    a=1/8.45; mu=MZ
    for thr,nf in ((4.18e9,5),(1.27e9,4),(0.0,3)):
        lo=max(E,thr)
        if lo<mu: a=1/(1/a+(11-2*nf/3)/(2*math.pi)*L0(lo/mu)); mu=lo
        if E>=thr: break
    return a
Es2=np.logspace(lg(4.5e8),lg(MZ),120); line([(X(lg(min(1.0,as_(E)))),Y(lg(E))) for E in Es2],ST)
line([(X(0),Y(lg(4.5e8))),(X(0),Y(lg(1.55e8)))],ST)
o.append(f'<circle cx="{X(0):.1f}" cy="{Y(lg(1.55e8)):.1f}" r="7" fill="{ST}" stroke="#fff" stroke-width="1.5"/>')
# weak below W: effective strength falls as (E/80 GeV)^2 ... draw as dashed fading
Es3=np.logspace(lg(1e3),lg(80.4e9),80); line([(X(lg(1/29.6)+2*lg(E/80.4e9)),Y(lg(E))) for E in Es3 if lg(1/29.6)+2*lg(E/80.4e9)>-45],WK,2.4,'stroke-dasharray="6 5"')
# gravity
Es4=np.logspace(lg(1e-6),LG_,200); line([(X(2*lg(E/1.2209e28)),Y(lg(E))) for E in Es4 if 2*lg(E/1.2209e28)>-45],GV,3)
IX0,IX1,IY0,IY1=GX0+70,GX0+500,Y(lg(1e-27)),Y(lg(1e-7))
o.append(f'<rect x="{IX0-50}" y="{IY1-46}" width="{IX1-IX0+90}" height="{IY0-IY1+100}" rx="8" fill="#f8fafc" stroke="#9aa7b2"/>')
txt(IX0-36,IY1-22,"Zoom: the three strengths from the Z up to where they meet",14,"#1f2a33",700,halo=False)
ix=lambda v: IX0+(IX1-IX0)*v/65; iy=lambda E: IY0-(IY0-IY1)*(lg(E)-lg(MZ))/(lg(4.62e22)-lg(MZ))
for v in (0,20,40,60): o.append(f'<line x1="{ix(v):.1f}" y1="{IY0}" x2="{ix(v):.1f}" y2="{IY1}" stroke="#e3e8ee"/>'); txt(ix(v),IY0+18,("1/%d"%v) if v else "full",12,"#56616b",a="middle",halo=False)
txt((IX0+IX1)/2,IY0+38,"strength (weaker to the right)",12.5,"#3b4650",a="middle",halo=False)
for i,col in ((0,EL),(1,WK),(2,ST)): line([(ix(gauge(E)[i]),iy(E)) for E in Es],col,3)
o.append(f'<line x1="{IX0}" y1="{iy(2.41e17):.1f}" x2="{IX1}" y2="{iy(2.41e17):.1f}" stroke="{HI}" stroke-dasharray="3 4"/>'); txt(IX1+4,iy(2.41e17)+4,"MMG",12.5,HI,700,halo=False)
txt(IX1+4,IY1+6,"meet",12.5,HI,700,halo=False); txt(IX1+4,IY0,"Z",12.5,"#56616b",700,halo=False)
txt(ix(8.45)+8,IY0-8,"strong",13,ST,700,halo=False); txt(ix(29.6)+8,IY0-8,"weak",13,WK,700,halo=False); txt(ix(59)-6,IY0-8,"electric side",13,EL,700,"end",halo=False)
# meeting marker
ym=Y(lg(4.62e22)); xm=X(-lg(39.0)); o.append(f'<circle cx="{xm:.1f}" cy="{ym:.1f}" r="9" fill="none" stroke="{HI}" stroke-width="2.5"/>')
for E in (2.41e17,4.62e22,80.4e9,1.55e8,3.73e3): o.append(f'<line x1="{GX0}" y1="{Y(lg(E)):.1f}" x2="{GX1}" y2="{Y(lg(E)):.1f}" stroke="#9aa7b2" stroke-width="1" stroke-dasharray="2 5"/>')
RX=GX1+16
def rl(l,t,s,col,dy=0):
    txt(RX,Y(l)+dy,t,14.5,col,700); 
    for i,q in enumerate(s): txt(RX,Y(l)+dy+18+i*16,q,12.5,"#3b4650")
rl(LG_,"gravity reaches full strength",["only at the grain"],GV,dy=10)
rl(lg(4.62e22),"the three meet (1/39)",["4.6 × 10¹³ GeV"],HI,dy=-2)
rl(lg(2.41e17),"the weak line bends",["at the MMG, 2.4 × 10⁸ GeV"],HI,dy=4)
rl(lg(80.4e9),"below the W the weak",["force fades fast: short reach"],WK,dy=4)
rl(lg(1.55e8),"strong reaches full strength",["and locks quarks into protons"],ST,dy=4)
rl(lg(3.73e3),"electric stays at 1/137",["all the way down: endless reach"],EL,dy=4)
txt(X(-36),Y(lg(1e3))-10,"gravity",15,GV,700,"middle"); txt(X(-13),Y(lg(3e6)),"weak (below the W)",14,WK,700,"middle")
txt(X(-2.1)-12,Y(lg(1e-10)),"electric",15,EL,700,"end"); txt(X(-0.9)-70,Y(lg(3e12)),"strong",15,ST,700,"end")
# legend panel
PY=Y0+RW*YS+60; o.append(f'<rect x="40" y="{PY}" width="{W-80}" height="{H-PY-26}" rx="10" fill="#fff" stroke="#c9d3dc"/>')
txt(64,PY+32,"Reading it",19,w=700,halo=False)
rows=[("Along a spiral","Nothing but the link that gives weight (the Higgs link). A rung's height is its grip. Strong, electric and weak treat all three rungs of a spiral alike, so they cannot tell the rungs apart."),
("Between spirals","The weak force. It turns an up-kind into a down-kind and an electron-kind into its neutrino. The small tilt between the up and down spirals is the quark mixing."),
("Around the first rungs","The strong and electric pulls, the two that bind. They make the proton ring and the atom ring. Only first-rung things last long enough to be ringed."),
("The cone itself","Gravity. It is not a mark on the cone; it is the slope from surface to grain. It is weaker than the others by 34 to 45 tens at the heights where particles sit, and equal only at the grain.")]
for i,(a,b) in enumerate(rows):
    y=PY+62+i*27; txt(64,y,a,14.5,"#1f2a33",700,halo=False); txt(270,y,b,14,"#3b4650",halo=False)
print(PY+62+4*27+20)
open("forces.svg","w").write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Inter, Helvetica, Arial, sans-serif"><rect width="{W}" height="{H}" fill="{BG}"/>'+"".join(o)+'</svg>')
