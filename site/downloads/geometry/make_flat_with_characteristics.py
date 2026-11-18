import math
F_='font-family="Inter, Helvetica, Arial, sans-serif"'
W,H=1150,1640; CX,CY,R=575,650,525
o=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" {F_}><rect width="{W}" height="{H}" fill="#f4f7fa"/>']
def txt(x,y,s,sz=15,c="#1f2a33",w=400,a="start",halo=True):
    h=' paint-order="stroke" stroke="#eaf0f5" stroke-width="4"' if halo else ''
    o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{sz}" fill="{c}" font-weight="{w}" text-anchor="{a}"{h}>{s}</text>')
txt(30,44,"The flat image with each particle's characteristics",27,w=700,halo=False)
txt(30,72,"Inside each circle: electric charge. Arc on top: spin. Small dots: colour. Beside it: energy at rest, grip on the surface, and how long it lasts.",15,"#3b4650",halo=False)
o.append(f'<circle cx="{CX}" cy="{CY}" r="{R}" fill="#e6eef5" stroke="#9fb3c4" stroke-width="2"/>')
def pol(r,deg): a=math.radians(deg); return CX+r*math.cos(a),CY-r*math.sin(a)
x1,y1=pol(R,180); x2,y2=pol(R,240)
o.append(f'<path d="M{CX},{CY} L{x1:.1f},{y1:.1f} A{R},{R} 0 0 0 {x2:.1f},{y2:.1f} Z" fill="#F5C518" fill-opacity="0.15"/>')
for d in (210,330,90):
    x,y=pol(R-6,d); o.append(f'<line x1="{CX}" y1="{CY}" x2="{x:.1f}" y2="{y:.1f}" stroke="#9fb3c4" stroke-width="1.2" stroke-dasharray="3 5"/>')
NU,LE,DN,UP,GO="#7a5aa6","#2f6db5","#2e8b57","#d1495b","#8a6d00"
def spin_arc(x,y,r,s,col="#1f2a33"):
    if s==0: return
    rr=r+6
    if s==0.5:
        o.append(f'<path d="M{x-rr:.1f},{y:.1f} A{rr},{rr} 0 0 1 {x+rr:.1f},{y:.1f}" fill="none" stroke="{col}" stroke-width="2.2"/>')
        o.append(f'<polygon points="{x+rr-5:.1f},{y-7:.1f} {x+rr+5:.1f},{y-7:.1f} {x+rr:.1f},{y+3:.1f}" fill="{col}"/>')
    else:
        o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rr}" fill="none" stroke="{col}" stroke-width="2.2"/>')
        o.append(f'<polygon points="{x+rr-5:.1f},{y-4:.1f} {x+rr+5:.1f},{y-4:.1f} {x+rr:.1f},{y+6:.1f}" fill="{col}"/>')
def part(x,y,r,col,chg,s,colour=False,lasts=True):
    o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{col}" fill-opacity="{0.95 if lasts else 0.30}" stroke="{col}" stroke-width="2"/>')
    o.append(f'<text x="{x:.1f}" y="{y+5.5:.1f}" font-size="{15 if len(chg)<3 else 13}" font-weight="700" fill="{"#fff" if lasts else "#1f2a33"}" text-anchor="middle">{chg}</text>')
    spin_arc(x,y,r,s)
    if colour:
        for k,cc in enumerate(("#e74c3c","#27ae60","#2e86de")): o.append(f'<circle cx="{x-9+k*9:.1f}" cy="{y+r+9:.1f}" r="3.4" fill="{cc}"/>')
RR={"nu":175,"le":275,"dn":375,"up":472}
P=[ # ray, slot, colour, name, charge, energy, grip, lasts-text, lasting?, has colour
 (210,"nu",NU,"neutrino 1","0","about 0.4 meV","","lasts",True,False),
 (210,"le",LE,"electron","−1","0.511 MeV","grip 0.0000029","lasts",True,False),
 (210,"dn",DN,"down","−⅓","4.7 MeV","grip 0.000027","lasts inside a skin",True,True),
 (210,"up",UP,"up","+⅔","2.2 MeV","grip 0.000012","lasts inside a skin",True,True),
 (330,"nu",NU,"neutrino 2","0","8.6 meV","","lasts",True,False),
 (330,"le",LE,"muon","−1","105.7 MeV","grip 0.00061","2.2 × 10⁻⁶ s",False,False),
 (330,"dn",DN,"strange","−⅓","93.5 MeV","grip 0.00054","about 10⁻⁸ s",False,True),
 (330,"up",UP,"charm","+⅔","1.27 GeV","grip 0.0073","about 10⁻¹² s",False,True),
 (90,"nu",NU,"neutrino 3","0","50 meV","","lasts",True,False),
 (90,"le",LE,"tau","−1","1.777 GeV","grip 0.010","2.9 × 10⁻¹³ s",False,False),
 (90,"dn",DN,"bottom","−⅓","4.18 GeV","grip 0.024","1.5 × 10⁻¹² s",False,True),
 (90,"up",UP,"top","+⅔","172.6 GeV","grip 0.99","5 × 10⁻²⁵ s",False,True)]
for deg,k,col,name,chg,en,grip,last,ok,cl in P:
    x,y=pol(RR[k],deg); part(x,y,21,col,chg,0.5,cl,ok)
    l2=(grip+" · " if grip else "")+last
    if deg==90: tx,ty,a=x+40,y-4,"start"
    elif deg==210: tx,ty,a=x+30,y+34,"start"
    elif k=="nu": tx,ty,a=x+34,y-30,"start"
    else: tx,ty,a=x-30,y+40,"end"
    txt(tx,ty,name+"  "+en,14.5,col,700,a); txt(tx,ty+17,l2,13,"#3b4650",400,a)
# neutrino triangle
T=[pol(RR["nu"],d) for d in (210,330,90)]
o.insert(8,f'<polygon points="{" ".join(f"{a:.1f},{b:.1f}" for a,b in T)}" fill="{NU}" fill-opacity="0.06" stroke="{NU}" stroke-width="1.1" stroke-dasharray="5 4"/>')
# centre carriers
C=[(0,0,"#e0a800","0",1,"light"),(88,52,"#715800","0",1,"gluon"),(78,155,"#715800","±1",1,"W"),(50,262,"#715800","0",1,"Z")]
for r,deg,col,chg,s,name in C:
    x,y=pol(r,deg); part(x,y,17,col,chg,s,name=="gluon",name in("light","gluon"))
    txt(x-30 if name=="Z" else x,y+5 if name=="Z" else y-30,name,13.5,"#715800",700,"end" if name=="Z" else "middle")
x,y=pol(118,330-0)  # Higgs placed off-ray
x,y=pol(76,350); part(x,y,17,"#5b6b78","0",0,False,False); txt(x+26,y+5,"Higgs",13.5,"#3b4650",700)
# shelled bubbles
for (r,deg,name,chg,en,last) in ((330,188,"proton","+1","938.3 MeV","lasts"),(455,186,"neutron","0","939.6 MeV","lasts in a nucleus; 15 min alone")):
    x,y=pol(r,deg); o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="30" fill="#b08900" fill-opacity="0.22" stroke="#715800" stroke-width="3.5"/>')
    o.append(f'<text x="{x:.1f}" y="{y+6:.1f}" font-size="16" font-weight="700" fill="#1f2a33" text-anchor="middle">{chg}</text>'); spin_arc(x,y,30,0.5)
    if name=="proton": txt(x,y-62,name+"  "+en,14.5,"#715800",700,"middle"); txt(x,y-45,last,13,"#3b4650",400,"middle")
    else: txt(x-34,y+58,name+"  "+en,14.5,"#715800",700,"start"); txt(x-34,y+75,"lasts in a nucleus;",13,"#3b4650",400,"start"); txt(x-34,y+91,"15 min alone",13,"#3b4650",400,"start")
txt(CX,CY-R-12,"THIRD RUNG",14,"#5b6b78",700,"middle"); txt(52,CY+R*0.5+110,"FIRST RUNG: what lasts",15,"#1f2a33",700); txt(W-52,CY+R*0.5+110,"SECOND RUNG",14,"#5b6b78",700,"end")
# ---------- key + table ----------
y0=1215
o.append(f'<rect x="20" y="{y0}" width="{W-40}" height="{H-y0-20}" fill="#ffffff" rx="10" stroke="#d5dde4"/>')
txt(40,y0+34,"How to read it",20,w=700,halo=False)
kx,ky=70,y0+78
part(kx,ky,17,LE,"−1",0.5); txt(kx+36,ky-2,"Half-turn arc: spin ½",14.5,w=700,halo=False); txt(kx+36,ky+16,"every particle of matter",13.5,"#3b4650",halo=False)
kx=330; part(kx,ky,17,"#715800","0",1); txt(kx+36,ky-2,"Full ring: spin 1",14.5,w=700,halo=False); txt(kx+36,ky+16,"every force carrier",13.5,"#3b4650",halo=False)
kx=560; part(kx,ky,17,"#5b6b78","0",0,False,False); txt(kx+36,ky-2,"No arc: spin 0",14.5,w=700,halo=False); txt(kx+36,ky+16,"the Higgs only",13.5,"#3b4650",halo=False)
kx=770; part(kx,ky,17,UP,"+⅔",0.5,True); txt(kx+36,ky-2,"Three dots: carries colour",14.5,w=700,halo=False); txt(kx+36,ky+16,"feels the strong pull",13.5,"#3b4650",halo=False)
txt(40,ky+62,"Solid circle = lasts. Pale circle = breaks down. Number inside = electric charge. Every particle of matter also has an opposite with the opposite charge.",14.2,"#1f2a33",halo=False)
ty=ky+104
txt(40,ty,"The carriers on the centre",16.5,w=700,halo=False)
rows=[("light (photon)","spin 1","charge 0","no colour","no energy at rest, no grip","lasts; carries the electric pull"),
      ("gluon (8 kinds)","spin 1","charge 0","colour","no energy at rest, no grip","never met alone; carries the strong pull"),
      ("W","spin 1","charge ±1","no colour","80.4 GeV","3 × 10⁻²⁵ s; changes one kind into another"),
      ("Z","spin 1","charge 0","no colour","91.2 GeV","3 × 10⁻²⁵ s"),
      ("Higgs","spin 0","charge 0","no colour","125.2 GeV","1.6 × 10⁻²² s; the surface's own ripple")]
X=[40,200,290,400,510,760]
for i,r in enumerate(rows):
    for j,s in enumerate(r): txt(X[j],ty+28+i*25,s,14,"#1f2a33",700 if j==0 else 400,halo=False)
txt(40,ty+28+5*25+16,"Spin, charge and colour are exact and the same on every rung. Only the energy (the grip) and the lifetime change from rung to rung.",14.8,"#1f2a33",700,halo=False)
o.append('</svg>'); open("flat2.svg","w").write("\n".join(o))
