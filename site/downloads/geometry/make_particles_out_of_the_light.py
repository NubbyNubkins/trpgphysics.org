import math
W=2300; H=2190; o=[]; BG="#f4f7fa"
def txt(x,y,s,sz=15,c="#1f2a33",w=400,a="start",halo=True):
    h=f' paint-order="stroke" stroke="{BG}" stroke-width="4"' if halo else ''
    o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{sz}" fill="{c}" font-weight="{w}" text-anchor="{a}"{h}>{s}</text>')
NU,LE,DN,UP,GO,GR,HI,RED,SF,LT="#7a5aa6","#2f6db5","#2e8b57","#d1495b","#715800","#56616b","#444b9a","#c0392b","#1b6f8a","#c98a00"
lg=math.log10
GRAIN,H0,T0=8.633e27,1.438e-33,2.349e-4; LS,LG_=lg(H0),lg(GRAIN)
CX=560; ZT=1560; YS=0.22; RW=470; RT=40; Y0=230+ZT
def zr(l): u=(l-LS)/(LG_-LS); return u*ZT, RT+(1-u)*(RW-RT)
def P(l,deg): z,r=zr(l); a=math.radians(deg); return CX+r*math.cos(a), Y0-z-r*math.sin(a)*YS
def ell(l,col,sw,dash="",fill="none",op=1): z,r=zr(l); o.append(f'<ellipse cx="{CX}" cy="{Y0-z:.1f}" rx="{r:.1f}" ry="{r*YS:.1f}" fill="{fill}" stroke="{col}" stroke-width="{sw}" stroke-opacity="{op}" {dash}/>')
txt(40,52,"Particles created out of the spiral: where, when and why each one drops out of the light",30,w=700,halo=False)
txt(40,84,"Before the light, the universe was cold and dark. The light began in the yellow band and stretches down the cone as space grows. Each time it falls past a particle's level, it can no longer make it,",16,"#3b4650",halo=False)
txt(40,108,"so that kind drops out of the light and what is left of it stays behind as a held thing. The time on each line is the age of the universe when the light passed that level.",16,"#3b4650",halo=False)
txt(40,132,"Light keeps stretching. Held things stop stretching where they dropped out, which is why they stay at their level while the light goes on down.",16,"#3b4650",halo=False)
o.append(f'<polygon points="{CX-RW},{Y0} {CX-RT},{Y0-ZT} {CX+RT},{Y0-ZT} {CX+RW},{Y0}" fill="#eef0f4"/>')
ell(LS,SF,3.5,fill="#d6e8ee")
o.append(f'<line x1="{CX-RW}" y1="{Y0}" x2="{CX-RT}" y2="{Y0-ZT}" stroke="#9aa7b2" stroke-width="1.3"/><line x1="{CX+RW}" y1="{Y0}" x2="{CX+RT}" y2="{Y0-ZT}" stroke="#9aa7b2" stroke-width="1.3"/>')
# the light: one spiral, six turns grain->surface; solid down to today, dashed below (future)
LT0=lg(T0)
def helix(l0,l1,dash,op):
    pts=[]; n=900
    for i in range(n+1):
        l=l0+(l1-l0)*i/n; pts.append(P(l,90+2160*(LG_-l)/(LG_-LS)))
    o.append(f'<polyline points="{" ".join(f"{x:.1f},{y:.1f}" for x,y in pts)}" fill="none" stroke="{LT}" stroke-width="3.2" stroke-opacity="{op}" stroke-linecap="round" {dash}/>')
LA,LB=lg(1e15),lg(1.3e17)
za,ra=zr(LA); zb,rb=zr(LB)
o.append(f'<polygon points="{CX-ra:.1f},{Y0-za:.1f} {CX-rb:.1f},{Y0-zb:.1f} {CX+rb:.1f},{Y0-zb:.1f} {CX+ra:.1f},{Y0-za:.1f}" fill="#F5C518" fill-opacity="0.35"/>')
helix(LG_,(LA+LB)/2,'stroke-dasharray="2 9"',0.35); helix((LA+LB)/2,LT0,"",0.95); helix(LT0,LS,'stroke-dasharray="3 7"',0.6)
txt(CX-rb-14,Y0-(za+zb)/2+5,"the light began in this band",14,GO,700,"end")
ell(LG_,GO,3.5,fill="#f1eee4")
# events: (label, energy eV, age, what happened, colour)
EV=[("THE GRAIN",GRAIN,"never reached","the smallest possible wave; the light did not begin this high",GO),
("D, coloured far-side parts",4.62e22,"never in the light","above where the light began",HI),
("hidden partner 3",6.30e19,"never in the light","handed over cold at the start; its uneven fall leaves more matter than antimatter",HI),
("hidden partner 2",1.29e18,"never in the light","handed over cold at the start; the light began at least ten times below here",HI),
("MMG, far-side weak pairs",2.41e17,"never in the light","just above where the light began; any made as the light built up decay at once",HI),
("hidden partner 1",7.25e14,"4.5 × 10⁻¹⁹ s","the only far-side thing the light made; the dark matter wave begins to swing near here",HI),
("top",172.6e9,"8 × 10⁻¹² s","heaviest front-side particle goes (stopped by temperature)",UP),
("the surface's resistance switches on",160e9,"1 × 10⁻¹¹ s","TEMPERATURE ALONE: the Higgs settles; every grip now gives a weight",RED),
("Higgs, Z, W",91.2e9,"2 to 4 × 10⁻¹¹ s","the weak carriers go; the weak force becomes short-reach",GR),
("bottom, tau, charm",1.78e9,"0.15 to 1.7 × 10⁻⁷ s","the heavy rungs go",DN),
("quarks lock into protons and neutrons",155e6,"24 millionths of a second","TEMPERATURE ALONE: held by the strong pull; 99% of the weight of everything you touch",RED),
("muon, strange",1.0e8,"58 to 84 millionths of a second","the middle rungs go",LE),
("neutrinos let go",1.0e6,"0.7 s","SPEED AGAINST GROWTH: space grows faster than they can meet anything",NU),
("electron",0.511e6,"3 s","electron and anti-electron pairs turn to light; the small leftover stays",LE),
("first nuclei",7e4,"4 minutes","COUNT OF LIGHT: binding is 2.2 MeV, but 1.6 billion bits of light per bit of matter delay it to 0.07 MeV",GR),
("majoron level",8.8e3,"5 hours","dark matter candidate (made earlier, cold; not made here)",GR),
("atoms form, light set free",0.26,"363 thousand years","COUNT OF LIGHT: binding is 13.6 eV, delayed to 0.3 eV by the same 1.6 billion; this is the glow we measure",RED),
("neutrino 3 slows",0.0223,"18 million years","it stops behaving like light and becomes a held thing",NU),
("neutrino 2 slows",0.00384,"259 million years","same",NU),
("TODAY: the light is here",T0,"13.8 billion years","52% of the way down; one wave is about 2 mm",LT),
("neutrino 1 slows",1.6e-4,"19.7 billion years (not yet)","still 44% movement today",NU),
("THE SURFACE",H0,"never quite","the stretch has a floor here (Gibbons and Hawking)",SF)]
n=len(EV); XL=CX+RW+150; ytop=Y0-ZT-10; ybot=Y0+18
for i,(nm,E,age,what,col) in enumerate(EV):
    l=lg(E); z,r=zr(l); yl=ytop+(ybot-ytop)*i/(n-1); yr_=Y0-z
    if 0<i<n-1: ell(l,col,1.6,'stroke-dasharray="2 4"',op=0.8)
    x0=CX+r
    o.append(f'<polyline points="{x0+4:.1f},{yr_:.1f} {XL-60:.1f},{yl-6:.1f} {XL-10:.1f},{yl-6:.1f}" fill="none" stroke="{col}" stroke-width="1.3" stroke-opacity="0.75"/>')
    o.append(f'<circle cx="{x0:.1f}" cy="{yr_:.1f}" r="4.5" fill="{col}" stroke="#fff" stroke-width="1.2"/>')
    big = i in (0,n-1,19)
    txt(XL,yl-10,nm,17 if big else 16,col,700)
    txt(XL+370,yl-10,age,16,"#1f2a33",700)
    txt(XL,yl+10,what,13.5,"#3b4650")
txt(XL,ytop-52,"WHAT DROPS OUT",13,"#6b7884",700); txt(XL+370,ytop-52,"AGE OF THE UNIVERSE",13,"#6b7884",700)
# today marker ring
z,r=zr(LT0); o.append(f'<ellipse cx="{CX}" cy="{Y0-z:.1f}" rx="{r:.1f}" ry="{r*YS:.1f}" fill="none" stroke="{LT}" stroke-width="3"/>')
txt(CX,Y0-z+r*YS+24,"the light today",15,LT,700,"middle")
txt(CX,Y0-ZT-RT*YS-16,"the grain",14,GO,700,"middle")
txt(CX,Y0+RW*YS+34,"the surface: one wave as long as everything we can see",15,SF,700,"middle")
# panel: stretch vs held
px,py,pw=1560,1700,700
PX=40; PY=Y0+RW*YS+70
o.append(f'<rect x="{PX}" y="{PY}" width="{W-80}" height="{H-PY-30}" rx="10" fill="#fff" stroke="#c9d3dc"/>')
txt(PX+24,PY+34,"What stretches and what does not, today",20,w=700,halo=False)
txt(PX+24,PY+58,"Yellow = free movement, which stretches with space. Dark = held energy (equilibrium), which does not.",14.5,"#3b4650",halo=False)
B=[("light",100,"100% movement: all of it stretches"),("neutrino 1",43.7,"44% movement: still half light-like"),("neutrino 2",0.19,"0.19% movement"),("neutrino 3",0.0056,"0.006% movement"),("electron",6.9e-8,"0.00000007% movement: held"),("proton",3.8e-11,"0.00000000004% movement: held")]
bw=300
for i,(nm,pc,s) in enumerate(B):
    col=i%3; row=i//3; x=PX+24+col*735; y=PY+84+row*44
    txt(x,y+17,nm,15.5,"#1f2a33",700,halo=False)
    o.append(f'<rect x="{x+100}" y="{y}" width="{bw}" height="24" rx="4" fill="#33414d"/>')
    wv=max(bw*pc/100,0 if pc<0.05 else 2)
    if wv>0: o.append(f'<rect x="{x+100}" y="{y}" width="{wv:.1f}" height="24" rx="4" fill="#F5C518"/>')
    txt(x+100+bw+12,y+17,s,13.5,"#3b4650",halo=False)
print(PY+84+2*44+30)
open("born.svg","w").write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Inter, Helvetica, Arial, sans-serif"><rect width="{W}" height="{H}" fill="{BG}"/>'+"".join(o)+'</svg>')
