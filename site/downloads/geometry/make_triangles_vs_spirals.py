import math,cmath
F_='font-family="Inter, Helvetica, Arial, sans-serif"'
W,H=2340,5400
o=[]
BG="#f4f7fa"
def txt(x,y,s,sz=15,c="#1f2a33",w=400,a="start",halo=True,hc=None):
    if hc is None: hc=BG
    h=f' paint-order="stroke" stroke="{hc}" stroke-width="4"' if halo else ''
    o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{sz}" fill="{c}" font-weight="{w}" text-anchor="{a}"{h}>{s}</text>')
NU,LE,DN,UP,GO,GR,HI,RED="#7a5aa6","#2f6db5","#2e8b57","#d1495b","#715800","#56616b","#444b9a","#c0392b"
nu=[3.6e-4,(3.6e-4**2+7.42e-5)**0.5,(3.6e-4**2+2.51e-3)**0.5]
le=[0.51099895e6,105.658e6,1776.93e6]; dn=[4.70e6,93.5e6,4183e6]; up=[2.16e6,1273e6,172570e6]
N=[le[i]**2/nu[i] for i in range(3)]
MMG,DD,GRAIN,MAJ=2.41e17,4.62e22,8.63e27,8.8e3
ANG=[210,330,450]
ZFRONT=[(-3.7,-1.1,30,108),(5.5,11.5,138,424),(14.5,20.0,452,522)]
ZFAR=[(-3.7,-1.1,26,72),(5.5,11.5,104,300),(14.5,28.0,330,522)]
Z=list(ZFRONT)
def rad(eV):
    l=math.log10(eV)
    for a,b,r1,r2 in Z:
        if a-0.01<=l<=b+0.01: return r1+(l-a)/(b-a)*(r2-r1)
    raise ValueError(l)
ZS,YS=0.80,0.30
PROJ=["flat"]
def P(cx,cy,r,deg):
    a=math.radians(deg)
    if PROJ[0]=="flat": return cx+r*math.cos(a),cy-r*math.sin(a)
    return cx+r*math.cos(a), cy-r*ZS-r*math.sin(a)*YS
def ring(cx,cy,r,col,sw,dash="",fill="none",op=1):
    if PROJ[0]=="flat": o.append(f'<circle cx="{cx}" cy="{cy}" r="{r:.1f}" fill="{fill}" stroke="{col}" stroke-width="{sw}" stroke-opacity="{op}" {dash}/>')
    else: o.append(f'<ellipse cx="{cx}" cy="{cy-r*ZS:.1f}" rx="{r:.1f}" ry="{r*YS:.1f}" fill="{fill}" stroke="{col}" stroke-width="{sw}" stroke-opacity="{op}" {dash}/>')
LAD=[("neutrinos",nu,NU,2.2),("leptons",le,LE,2.2),("down quarks",dn,DN,3),("up quarks",up,UP,3),("hidden partners",N,HI,2.2)]
def fac(f): return ("×%.0f"%f) if f>=20 else "×%.1f"%f
TL={"neutrinos":(.5,.5),"leptons":(.36,.30),"down quarks":(.68,.70),"up quarks":(.52,.50),"hidden partners":(.5,.5)}
NAMES={"neutrinos":[("neutrino 1","0.36 meV (assumed)"),("neutrino 2","8.6 meV"),("neutrino 3","50 meV")],
 "leptons":[("electron","0.511 MeV"),("muon","105.7 MeV"),("tau","1.777 GeV")],
 "down quarks":[("down","4.7 MeV"),("strange","93.5 MeV"),("bottom","4.18 GeV")],
 "up quarks":[("up","2.2 MeV"),("charm","1.27 GeV"),("top","172.6 GeV")],
 "hidden partners":[("hidden 1","7 × 10⁵ GeV"),("hidden 2","1.3 × 10⁹ GeV"),("hidden 3","6.3 × 10¹⁰ GeV")]}
LBF={"neutrinos":[(-10,-8,"end"),(12,16,"start"),(-10,-6,"end")],
 "leptons":[(-12,-10,"end"),(-18,26,"end"),(-14,5,"end")],
 "down quarks":[(-12,20,"end"),(18,-18,"start"),(-14,5,"end")],
 "up quarks":[(-12,-10,"end"),(16,2,"start"),(-14,5,"end")],
 "hidden partners":[(-12,-10,"end"),(14,18,"start"),(-14,6,"end")]}
LBC={"neutrinos":[(-12,6,"end"),(12,6,"start"),(-10,-8,"end")],
 "leptons":[(-14,16,"end"),(16,-14,"start"),(-14,8,"end")],
 "down quarks":[(-14,-12,"end"),(16,22,"start"),(14,-4,"start")],
 "up quarks":[(-14,6,"end"),(16,6,"start"),(-14,4,"end")],
 "hidden partners":[(-14,4,"end"),(14,4,"start"),(14,-6,"start")]}
def wave(m):
    sq=[math.sqrt(v) for v in m]; A=sum(sq)/3; c=[q/A-1 for q in sq]
    z=(2/3)*sum(c[i]*cmath.exp(-1j*math.radians(120*i)) for i in range(3)); return A,abs(z),math.degrees(cmath.phase(z))
def panel(cx,cy,proj,mode,far,title,sub):
    global BG,Z
    PROJ[0]=proj; Z=list(ZFAR if far else ZFRONT); BG="#ece9f3" if far else "#eef2e6"; cone=(proj=="cone")
    RB=536; edge="#444b9a" if far else "#b08900"
    top=cy-RB*(ZS+YS) if cone else cy-RB
    txt(cx-RB,top-40,title,22,w=700,halo=False); txt(cx-RB,top-16,sub,14.5,"#3b4650",halo=False)
    if cone:
        o.append(f'<polygon points="{cx},{cy} {cx-RB},{cy-RB*ZS:.1f} {cx+RB},{cy-RB*ZS:.1f}" fill="{BG}" stroke="none"/>')
        ring(cx,cy,RB,edge,1.5,fill=BG)
        o.append(f'<line x1="{cx}" y1="{cy}" x2="{cx-RB}" y2="{cy-RB*ZS:.1f}" stroke="#9aa7b2" stroke-width="1.2"/><line x1="{cx}" y1="{cy}" x2="{cx+RB}" y2="{cy-RB*ZS:.1f}" stroke="#9aa7b2" stroke-width="1.2"/>')
        o.append(f'<line x1="{cx-RB-26}" y1="{cy}" x2="{cx-RB-26}" y2="{cy-RB*ZS:.1f}" stroke="#3b4650" stroke-width="1.5"/><polygon points="{cx-RB-26},{cy-RB*ZS-10:.1f} {cx-RB-31},{cy-RB*ZS:.1f} {cx-RB-21},{cy-RB*ZS:.1f}" fill="#3b4650"/>')
        txt(cx-RB-20,cy-20,"height = energy",12.5,"#3b4650",600,"start",True,"#f4f7fa")
    else: ring(cx,cy,RB,edge,1.5,fill=BG)
    for (a_,b_,r1,r2),(a2,b2,r3,r4) in zip(Z[:-1],Z[1:]):
        if cone: ring(cx,cy,r2,"#9aa7b2",1,'stroke-dasharray="2 5"'); ring(cx,cy,r3,"#9aa7b2",1,'stroke-dasharray="2 5"')
        else:
            ring(cx,cy,(r2+r3)/2,"#ffffff",r3-r2-8); ring(cx,cy,(r2+r3)/2,"#9aa7b2",1,'stroke-dasharray="2 5"')
    if not cone:
        x,y=P(cx,cy,(Z[1][3]+Z[2][2])/2,165); txt(x,y+4,"empty stretch (drawn shortened)",11,"#6b7884",600,"middle",True,"#ffffff")
    for e_,lab in ((1e6,"1 MeV"),(1e7,""),(1e8,"100 MeV"),(1e9,"1 GeV"),(1e10,""),(1e11,"100 GeV")):
        r=rad(e_); ring(cx,cy,r,"#c9d3dc",0.9)
        if lab and not cone and not far: x,y=P(cx,cy,r,172); txt(x,y+4,lab,11,"#8a97a3",400,"middle")
    for d in (210,330,90):
        x,y=P(cx,cy,RB,d); o.append(f'<line x1="{cx}" y1="{cy}" x2="{x:.1f}" y2="{y:.1f}" stroke="#9fb3c4" stroke-width="1.1" stroke-dasharray="3 5"/>')
    rw1,rw2=rad(80.4e9),rad(125.2e9); ring(cx,cy,(rw1+rw2)/2,GO,max(rw2-rw1,4),op=0.22); ring(cx,cy,rad(1.7e9),GO,1.3,'stroke-dasharray="3 5"')
    la=180 if cone else 128; an="end" if cone else "middle"; lx=-8 if cone else 0; hc="#f4f7fa" if cone else None
    if not far:
        ring(cx,cy,rad(MMG),RED,1.7,'stroke-dasharray="5 4"')
        for e_,lab,col in ((MMG,"MMG",RED),((80.4e9*125.2e9)**.5,"W, Z, Higgs",GO),(1.7e9,"gluon-only bubbles",GO)):
            x,y=P(cx,cy,rad(e_),la); txt(x+lx,y+4,lab,12.5,col,700,an,True,hc)
    o.append(f'<circle cx="{cx}" cy="{cy}" r="6" fill="#e0a800" stroke="#1f2a33" stroke-width="1.5"/>')
    if cone: txt(cx,cy+22,"light + gluons (the point of the cone: zero)",12,"#8a6d00",700,"middle",True,"#f4f7fa")
    elif not far: txt(cx+10,cy-9,"light + gluons",11.5,"#8a6d00",700,"start")
    LB=LBC if cone else LBF
    if far:   # far-side furniture and the links to the front
        ring(cx,cy,rad(GRAIN),GO,2.4,'stroke-dasharray="6 4"'); ring(cx,cy,rad(DD),GR,2,'stroke-dasharray="5 4"'); ring(cx,cy,rad(MMG),RED,2,'stroke-dasharray="5 4"')
        rmj=(Z[0][3]+Z[1][2])/2; ring(cx,cy,rmj,NU,1.6,'stroke-dasharray="2 4"')
        for i in range(3):   # mirror line neutrino - lepton - hidden on each rung line
            x1,y1=P(cx,cy,rad(nu[i]),ANG[i]); x2,y2=P(cx,cy,rad(N[i]),ANG[i]); o.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{NU}" stroke-width="5" stroke-opacity="0.30" stroke-linecap="round"/>')
            for e_,col,sh in ((MMG,RED,"c"),(DD,GR,"s")):
                x,y=P(cx,cy,rad(e_),ANG[i])
                if sh=="c": o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7.5" fill="#fff" stroke="{col}" stroke-width="2.4"/>')
                else: o.append(f'<rect x="{x-7:.1f}" y="{y-7:.1f}" width="14" height="14" fill="#fff" stroke="{col}" stroke-width="2.4"/>')
        # second mirror: MMG - D - grain, drawn along a free direction
        am=150; x1,y1=P(cx,cy,rad(MMG),am); x2,y2=P(cx,cy,rad(GRAIN),am); xm,ym=P(cx,cy,rad(DD),am)
        o.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{GR}" stroke-width="5" stroke-opacity="0.45" stroke-linecap="round"/><circle cx="{xm:.1f}" cy="{ym:.1f}" r="5" fill="{GR}"/>')
        txt(xm+12,ym+22,"MMG × grain = D²",13.5,GR,700,"start"); txt(xm+12,ym+39,"(D is the mirror point)",12,"#3b4650",400,"start")
        # D presses on the down ladder (third rung)
        xa,ya=P(cx,cy,rad(DD),90); xb,yb=P(cx,cy,rad(dn[2])+12,90)
        o.append(f'<line x1="{xa+12:.1f}" y1="{ya+8:.1f}" x2="{xb+12:.1f}" y2="{yb:.1f}" stroke="{DN}" stroke-width="2.2" stroke-dasharray="6 4"/><polygon points="{xb+12:.1f},{yb+9:.1f} {xb+6:.1f},{yb-3:.1f} {xb+18:.1f},{yb-3:.1f}" fill="{DN}"/>')
        xt,yt=P(cx,cy,(rad(DD)+rad(N[2]))/2-6,90); 
        la2=[(GRAIN,"the grain 8.6 × 10¹⁸ GeV",GO),(DD,"D: coloured far-side parts, 4.6 × 10¹³ GeV (three, split unknown)",GR),(MMG,"MMG 2.4 × 10⁸ GeV (three, split unknown)",RED)]
        for e_,lab,col in la2:
            x,y=P(cx,cy,rad(e_),58 if not cone else 0); txt(x+(-12 if cone else 0),y-8,lab,13,col,700,"end" if cone else "middle")
        x,y=P(cx,cy,rmj,58 if not cone else 0); txt(x+(-8 if cone else 0),y+(-6 if cone else 4),"majoron 8.8 keV",12,NU,700,"end" if cone else "middle")
    for name,m,col,sw in LAD:
        hid=(name=="hidden partners"); dim=far and not hid; op=0.5 if dim else 0.9
        if far and hid: sw=3.4
        r=[rad(v) for v in m]; pts=[P(cx,cy,r[i],ANG[i]) for i in range(3)]
        dash=' stroke-dasharray="7 4"' if (hid and not far) else ''
        zone=[z_ for z_ in Z if z_[0]-0.01<=math.log10(m[0])<=z_[1]+0.01][0]
        d1,d2=r[1]-r[0],r[2]-r[1]
        def rq(t): return r[0]+d1*t+(d2-d1)*t*(t-1)/2
        if mode=="tri":
            o.append(f'<polygon points="{" ".join(f"{a_:.1f},{b_:.1f}" for a_,b_ in pts)}" fill="{col}" fill-opacity="{0.02 if dim else 0.05}" stroke="{col}" stroke-width="{sw}" stroke-opacity="{op}"{dash}/>')
        elif mode=="spi":
            path=[]
            for j in range(0,241,3):
                x,y=P(cx,cy,rq(j/120),210+j); path.append(f'{"M" if j==0 else "L"}{x:.1f},{y:.1f}')
            o.append(f'<path d="{" ".join(path)}" fill="none" stroke="{col}" stroke-width="{sw+0.6}" stroke-opacity="{op}" stroke-linecap="round"{dash}/>')
        else:
            A,e,th0=wave(m); path=[]
            for j in range(0,361):
                E=(A*(1+e*math.cos(math.radians(th0+j))))**2; l=max(math.log10(max(E,1e-300)),zone[0])
                rr=zone[2]+(l-zone[0])/(zone[1]-zone[0])*(zone[3]-zone[2]); x,y=P(cx,cy,rr,210+j); path.append(f'{"M" if j==0 else "L"}{x:.1f},{y:.1f}')
            o.append(f'<path d="{" ".join(path)}Z" fill="{col}" fill-opacity="{0.02 if dim else 0.04}" stroke="{col}" stroke-width="{sw+0.4}" stroke-opacity="{op}" stroke-linejoin="round"{dash}/>')
            if not dim:
                x,y=P(cx,cy,rad((A*(1+e))**2),210-th0); o.append(f'<rect x="{x-6:.1f}" y="{y-6:.1f}" width="12" height="12" transform="rotate(45 {x:.1f} {y:.1f})" fill="#fff" stroke="{col}" stroke-width="2.2"/>')
        if mode!="wav" and not dim:
            for s_ in range(2):
                t=TL[name][s_]
                if mode=="tri": x=pts[s_][0]+(pts[s_+1][0]-pts[s_][0])*t; y=pts[s_][1]+(pts[s_+1][1]-pts[s_][1])*t
                else: x,y=P(cx,cy,rq(s_+t),210+120*(s_+t))
                txt(x,y+5,fac(m[s_+1]/m[s_]),13.5,col,700,"middle")
        for i,(x,y) in enumerate(pts):
            if hid and not far: o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="8" fill="#fff" stroke="{col}" stroke-width="2.4" stroke-dasharray="4 3"/>')
            elif dim: o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{col}" fill-opacity="0.7"/>')
            else: o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{5.5 if name=="neutrinos" else (9 if hid else 7.5)}" fill="{col}" stroke="#fff" stroke-width="1.5"/>')
            if not dim:
                dx,dy,a_=LB[name][i]
                if far and hid: dx,dy,a_=[(-16,18,"end"),(16,20,"start"),(-16,-4,"end")][i]
                n1,n2=NAMES[name][i]; txt(x+dx,y+dy,n1+" "+n2,13 if hid and far else 12.5,col,700,a_)
    BG="#f4f7fa"
def head(y,t,s,c="#1f2a33"):
    txt(40,y,t,24,c,700,halo=False); txt(40,y+26,s,15.5,"#3b4650",halo=False)
txt(40,48,"The ladders as triangles, spirals and waves: front side and far side, flat and with height",28,w=700,halo=False)
txt(40,76,"Distance from the centre (or height) is energy; every equal step is ten times more. Same particles and same × steps throughout. Only the way they are joined changes.",15.5,"#3b4650",halo=False)
RB=536; PH=RB*(ZS+YS); y=110
head(y+20,"1. Front side, flat","The particles we have met, plus the hidden-partner ladder (dashed) for reference. Left: joined by straight lines. Right: joined by turning outward.")
c=y+20+26+70+RB
panel(585,c,"flat","tri",False,"Front side: triangles","Each ladder: three rungs, one third of a turn apart."); panel(1755,c,"flat","spi",False,"Front side: spirals","Each ladder: turn a third of a turn, step out, turn again, step out.")
y=c+RB+60
head(y,"2. Far side, flat, and how it links to the front","Far-side particles are bold. The front ladders are kept faint in the middle so the links show. The outer band is squeezed to fit the far side in.",HI)
c=y+26+70+RB
panel(585,c,"flat","tri",True,"Far side: triangles","Only the hidden partners have three known rungs. The MMGs and Ds are three each, but only their average is known."); panel(1755,c,"flat","spi",True,"Far side: spirals","The hidden-partner spiral, with the same links.")
y=c+RB+60
head(y,"3. Front side, with height","The flat chart lifted into a cone. Height is energy, exactly as distance from the centre was.")
c=y+26+70+PH
panel(585,c,"cone","tri",False,"Front side, with height: triangles","The same three rungs per ladder, now climbing."); panel(1755,c,"cone","spi",False,"Front side, with height: spirals","Each ladder turns and climbs.")
y=c+70
head(y,"4. Far side, with height","The far side is the upper part of the same cone. The front ladders sit low in it, faint.",HI)
c=y+26+70+PH
panel(585,c,"cone","tri",True,"Far side, with height: triangles","Thick purple lines: each neutrino, its lepton and its hidden partner lie on one line."); panel(1755,c,"cone","spi",True,"Far side, with height: spirals","The hidden-partner spiral climbing above the front ladders.")
y=c+70
head(y,"5. The connections calculated as waves, with height","Each ladder is ONE wave going once round the cone; its three particles are three points on it. Diamond = the crest. Where a wave drops to the floor of its band it is touching zero.")
c=y+26+70+PH
panel(585,c,"cone","wav",False,"Front side: waves","The crest sits just before the heaviest rung. The lightest rung sits just before the wave touches zero."); panel(1755,c,"cone","wav",True,"Far side: waves","The hidden-partner wave (bold) above the front waves (faint).")
BG="#ffffff"; ty=c+60
def fE(E):
    for u_,f_ in (("GeV",1e9),("MeV",1e6),("keV",1e3),("eV",1),("meV",1e-3)):
        if E>=f_:
            v=E/f_; return (("%.3g"%v) if v<1e4 else ("%.1f × 10"%(v/10**int(math.log10(v)))+str(int(math.log10(v))).translate(str.maketrans("0123456789","⁰¹²³⁴⁵⁶⁷⁸⁹"))))+" "+u_
    return "%.2g meV"%(E/1e-3)
tb=len(o); o.append("")
txt(40,ty+34,"How the far side and the front side interact",20,w=700,halo=False)
L1=["1. Mirror on every rung (exact in the theory): neutrino × hidden partner = lepton². The lepton is the mirror point; that is the thick purple line on each of the three directions.",
"2. A second mirror of the same kind (exact in the theory): MMG × grain = D². Check: (4.62 × 10¹³)² ÷ 8.63 × 10¹⁸ = 2.47 × 10⁸ GeV, against 2.41 × 10⁸ GeV used for the MMG, and 2.2 × 10⁸ GeV needed to close the force strengths.",
"3. D presses on the down-quark ladder (green dashed arrow). The size of the press, 34.5, is measured, not derived.",
"4. Noticed: the hidden ladder and the up-quark ladder span almost the same total distance (÷86,887 and ÷79,894). They are equal if the lightest neutrino is 0.33 meV.",
"5. Noticed: the MMG sits near the middle of hidden 1 and hidden 3 (2.1 × 10⁸ against 2.4 × 10⁸ GeV). The majoron comes from the hidden partners settling."]
for i,l in enumerate(L1): txt(40,ty+66+i*25,l,15,halo=False)
y2=ty+66+5*25+26
txt(40,y2,"How much each spiral tightens (second step compared with the first; 1.00 would be a perfectly even spiral)",19,w=700,halo=False)
rows=[("Hidden partners","×1785, then ×49","0.52",HI),("Leptons","×207, then ×16.8","0.53",LE),("Neutrinos","×24, then ×5.8","0.55 (if neutrino 1 is 0.36 meV)",NU),("Up quarks","×589, then ×136","0.77",UP),("Down quarks","×19.9, then ×44.7","1.27 (opens instead of tightening)",DN)]
for i,(a_,b_,c_,col) in enumerate(rows):
    x=40+(i%3)*760; yy=y2+34+(i//3)*28
    txt(x,yy,a_,15,col,700,halo=False); txt(x+150,yy,b_,15,halo=False); txt(x+330,yy,c_,15,"#1f2a33",700,halo=False)
wy=y2+34+2*28+22
txt(40,wy,"The waves, calculated",19,w=700,halo=False)
hdr=["Ladder","crest before heaviest rung (the twist)","lightest rung before the zero","wave height ÷ its average","small rise between the two zeros"]
XX=[40,230,640,960,1280]
for j,h_ in enumerate(hdr): txt(XX[j],wy+30,h_,14,"#3b4650",700,halo=False)
for i,(name,m,col,sw) in enumerate(LAD):
    A,e,th0=wave(m); zero=math.degrees(math.acos(-1/e)); yy=wy+56+i*25
    for j,v in enumerate([name.capitalize(),"%.2f°"%abs(th0+240-360),"%.2f°"%(zero-th0),"%.3f"%e,fE((A*(1-e))**2)]): txt(XX[j],yy,v,14.5,col if j==0 else "#1f2a33",700 if j==0 else 400,halo=False)
yy=wy+56+5*25+10
txt(40,yy,"Three points always fit one wave, so the wave describes a ladder; it does not yet predict it. The spiral and wave lines are drawn through the measured points. Hidden-partner values assume the lightest neutrino is 0.36 meV.",14.5,"#3b4650",halo=False)
H=int(yy+40)
o[tb]=f'<rect x="20" y="{ty}" width="{W-40}" height="{H-ty-20}" fill="#ffffff" rx="10" stroke="#d5dde4"/>'
out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" {F_}><rect width="{W}" height="{H}" fill="#f4f7fa"/>']+o+['</svg>']
open("spiral3.svg","w").write("\n".join(out)); print(H, (4.62e13)**2/8.63e18, (N[0]*N[2])**.5/1e9, N[2]/N[0], up[2]/up[0])
