from onepic import *
F_='font-family="Inter, Helvetica, Arial, sans-serif"'
W,HG,HL=1150,1180,470; H=HG+HL+70
VIEW=(-35,22); r,u,v=camera(VIEW)
pts=[CORN[k] for k in CORN]+[np.array([0,0,HTOP+1.5])]
xs=[p@r for p in pts]; ys=[p@u for p in pts]
SC=0.80*min(W/(max(xs)-min(xs)+4),HG/(max(ys)-min(ys)+3)); cx=(max(xs)+min(xs))/2; cy=(max(ys)+min(ys))/2
def P(p): p=np.asarray(p,float); return (W/2+SC*(p@r-cx), 70+HG/2-SC*(p@u-cy))
def ax(eV,dx=0): return np.array([dx,0,hgt(eV)])
def txt(x,y,s,sz=15,c="#1f2a33",w=400,a="start",halo=True):
    h=' paint-order="stroke" stroke="#f4f7fa" stroke-width="4"' if halo else ''
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{sz}" fill="{c}" font-weight="{w}" text-anchor="{a}"{h}>{s}</text>'
def line(a,b,st): (x1,y1),(x2,y2)=P(a),P(b); return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" {st}/>'
def frame(st_near,st_far,axis_st):
    o=[]
    for a,b in EDGES:
        far=(sum(a)+sum(b))>=3; o.append(line(CORN[a],CORN[b],st_far if far else st_near))
    o.append(line(np.zeros(3),np.array([0,0,HTOP]),axis_st)); return o
def head(title,sub):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" {F_}><rect width="{W}" height="{H}" fill="#f4f7fa"/>',
            txt(30,44,title,27,w=700,halo=False),txt(30,72,sub,15.5,"#3b4650",halo=False)]
COLF={n:c for n,c,_ in LAD}
# ================= IMAGE 1: THINGS =================
o=head("Image 1. Things: what can be experienced","Only what leaves a mark is drawn. The faint outline is kept just so you can see where each thing sits.")
g='stroke="#c9d3dc" stroke-width="1" stroke-dasharray="2 5"'
o+=frame(g,g,g)
def bubble(p,rad,col,kind,lab,dx=0,dy=0,anchor="start"):
    x,y=P(p); s=''
    if kind=="A": s=f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rad}" fill="{col}" stroke="#1f2a33" stroke-width="2.2"/>'
    if kind=="B": s=f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rad}" fill="{col}" fill-opacity="0.85" stroke="#fff" stroke-width="1.5"/><circle cx="{x:.1f}" cy="{y:.1f}" r="{rad+5}" fill="none" stroke="{col}" stroke-width="1.5"/>'
    if kind=="C": s=f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rad}" fill="{col}" fill-opacity="0.30" stroke="{col}" stroke-width="1.6"/>'
    if kind=="D": s=f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{rad}" fill="none" stroke="#8a97a3" stroke-width="1.8" stroke-dasharray="4 4"/>'
    tx=x+(rad+7 if anchor=="start" else -rad-7)+dx
    return s+txt(tx,y+5+dy,lab,14.5,col if kind!="D" else "#6b7884",600,anchor)
GOLD="#b08900"
# lasting, direct
o.append(bubble(POS("e"),10,COLF["charged leptons"],"A","electron",anchor="end"))
o.append(bubble(POS("p"),15,GOLD,"A","proton, neutron, atoms (bubbles with a skin)",anchor="end"))
o.append(bubble(np.zeros(3),11,"#e0a800","A","light (photon)",dy=16))
# brief track
o.append(bubble(POS("μ"),9,COLF["charged leptons"],"B","muon",anchor="end"))
o.append(bubble(POS("τ"),8,COLF["charged leptons"],"B","tau",anchor="end",dy=4))
nuc=sum(ppos("neutrinos",k,m) for k,m in enumerate(NU))/3
o.append(bubble(nuc,9,COLF["neutrinos"],"B","neutrinos (rare clicks, always as a mix of three)"))
# only by what they leave
for lab,nm,an,dy in (("u","up","start",10),("d","down","start",-4),("s","strange","start",8),("c","charm","start",12),("b","bottom","start",-2),("t","top","start",0)):
    nmf="up-type" if lab in "uct" else "down-type"
    o.append(bubble(POS(lab),8,COLF[nmf],"C",nm,anchor=an,dy=dy))
o.append(bubble(ax(1.0e11),9,"#715800","C","W, Z, Higgs",anchor="end"))
o.append(bubble(np.array([1.6,0,0.9]),8,"#715800","C","gluons",dy=0))
# claimed, not yet experienced
for eV,lab,an in ((8.8e3,"majoron (dark matter candidate)","start"),(2.47e17,"far-side weak pairs","start"),(1e19,"hidden neutrino partners","end"),(4.62e22,"coloured far-side parts (D)","start"),(8.63e27,"the grain","start")):
    o.append(bubble(ax(eV),9,"#8a97a3","D",lab,anchor=an))
# pulls
bx=W-330
for i,(t1,t2) in enumerate((("gravity","felt as weight"),("dark energy","felt as space stretching"),("dark matter's pull","felt in how galaxies turn"))):
    y=150+i*62
    o.append(f'<path d="M{bx},{y} q10,-12 20,0 t20,0 t20,0" fill="none" stroke="#5b6b78" stroke-width="2.4"/><polygon points="{bx+60},{y-6} {bx+72},{y} {bx+60},{y+6}" fill="#5b6b78"/>')
    o.append(txt(bx+82,y-2,t1,15.5,"#1f2a33",700)); o.append(txt(bx+82,y+17,t2,13.5,"#3b4650"))
o.append(txt(bx,112,"Experienced as a pull, not as an object",15,"#1f2a33",700))
# legend
y0=70+HG+10
o.append(f'<rect x="20" y="{y0}" width="{W-40}" height="{HL-20}" fill="#ffffff" rx="10" stroke="#d5dde4"/>')
o.append(txt(40,y0+34,"How each thing is experienced",20,w=700,halo=False))
rows=[("A","#2f6db5","Directly, and it lasts","electron, proton, neutron, atoms, light. Antimatter too (positrons are used in hospital PET scans).",
       "The proton and the atom have a skin (the balance point). The electron shows no skin: it is felt as a single click and as the cloud round an atom."),
      ("B","#2f6db5","Directly, but only briefly or rarely","muon and tau leave a short track. Neutrinos give a rare click.","Neutrinos are only ever met as a mix of three, never one pure rung."),
      ("C","#2e8b57","Only by what it leaves behind","quarks, gluons, W, Z, Higgs, top.","Quarks and gluons are never met alone: they are met as the inside of the proton bubble, or as a spray. W, Z, Higgs and top break up before reaching anything."),
      ("D","#8a97a3","Claimed by the theory, not yet experienced","majoron, far-side weak pairs, hidden neutrino partners, coloured far-side parts, the grain.","These are meant as things, not concepts. Nobody has met one yet.")]
for i,(k,col,t1,t2,t3) in enumerate(rows):
    yy=y0+74+i*92; x=58
    if k=="A": o.append(f'<circle cx="{x}" cy="{yy}" r="11" fill="{col}" stroke="#1f2a33" stroke-width="2.2"/>')
    if k=="B": o.append(f'<circle cx="{x}" cy="{yy}" r="9" fill="{col}" stroke="#fff" stroke-width="1.5"/><circle cx="{x}" cy="{yy}" r="14" fill="none" stroke="{col}" stroke-width="1.5"/>')
    if k=="C": o.append(f'<circle cx="{x}" cy="{yy}" r="11" fill="{col}" fill-opacity="0.30" stroke="{col}" stroke-width="1.6"/>')
    if k=="D": o.append(f'<circle cx="{x}" cy="{yy}" r="11" fill="none" stroke="{col}" stroke-width="1.8" stroke-dasharray="4 4"/>')
    o.append(txt(88,yy-8,t1,16.5,w=700,halo=False)); o.append(txt(88,yy+13,t2,14.5,"#1f2a33",halo=False)); o.append(txt(88,yy+33,t3,13.5,"#3b4650",halo=False))
o.append('</svg>'); open("things.svg","w").write("\n".join(o))
# ================= IMAGE 2: CONCEPTS =================
o=head("Image 2. Concepts: the relationship chart","Nothing here can be met on its own. These are the relationships and rules that say how the things in Image 1 behave.")
o+=frame('stroke="#3b4650" stroke-width="2"','stroke="#8e3b8e" stroke-width="1.4" stroke-dasharray="6 4"','stroke="#715800" stroke-width="1.8" stroke-dasharray="3 5"')
T=[CORN[(0,0,0)],CORN[(1,1,0)],CORN[(1,0,1)],CORN[(0,1,1)]]; TM=[CORN[(1,1,1)],CORN[(1,0,0)],CORN[(0,1,0)],CORN[(0,0,1)]]
for TT,fill,st in ((TM,"#5aa0c8",'stroke="#2b7bb9" stroke-width="1.4" stroke-dasharray="7 4"'),(T,"#F5C518",'stroke="#b08900" stroke-width="1.4"')):
    for f in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]:
        Q=[P(TT[i]) for i in f]; o.append(f'<polygon points="{" ".join(f"{a:.1f},{b:.1f}" for a,b in Q)}" fill="{fill}" fill-opacity="0.06" {st}/>')
hc=6.0
for name,col,_ in LAD:
    rr=hc*math.tan(OPEN[name]*D2R)
    ring=[np.array([rr*math.cos(t),rr*math.sin(t),hc]) for t in np.linspace(0,2*math.pi,73)]
    o.append(f'<polyline points="{" ".join(f"{a:.1f},{b:.1f}" for a,b in map(P,ring))}" fill="none" stroke="{col}" stroke-width="{2.6 if name=="charged leptons" else 1.4}" stroke-opacity="0.8"/>')
    for t in np.linspace(0,2*math.pi,7)[:-1]:
        o.append(line(np.zeros(3),np.array([rr*math.cos(t),rr*math.sin(t),hc]),f'stroke="{col}" stroke-width="0.8" stroke-opacity="0.35"'))
for name,col,it in LAD:
    P3=[ppos(name,k,m) for k,(_,m) in enumerate(it)]
    o.append(f'<polyline points="{" ".join(f"{a:.1f},{b:.1f}" for a,b in map(P,P3))}" fill="none" stroke="{col}" stroke-width="2" stroke-opacity="0.8"/>')
    for p in P3:
        x,y=P(p); o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="#fff" stroke="{col}" stroke-width="1.6"/>')
# meeting X
xm,ym=P(ax(4.62e22))
for dxx,c in ((-70,"#c0392b"),(0,"#2b7bb9"),(70,"#2e8b57")):
    o.append(f'<line x1="{xm+dxx}" y1="{ym+150}" x2="{xm-dxx*0.5}" y2="{ym-75}" stroke="{c}" stroke-width="1.6" stroke-opacity="0.75"/>')
# centre
xc,yc=P(np.array([0,0,HTOP/2])); o.append(f'<circle cx="{xc:.1f}" cy="{yc:.1f}" r="6" fill="#16a085"/>')
lr=hc*math.tan(45*D2R)
def badge(p,n,dx=0,dy=0):
    x,y=P(p); x+=dx; y+=dy
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="13" fill="#1f2a33"/><text x="{x:.1f}" y="{y+5:.1f}" font-size="14" font-weight="700" fill="#fff" text-anchor="middle">{n}</text>'
up=[ppos("up-type",k,m) for k,(_,m) in enumerate(LAD[3][2])]; dn=[ppos("down-type",k,m) for k,(_,m) in enumerate(LAD[2][2])]
lp=[ppos("charged leptons",k,m) for k,(_,m) in enumerate(LAD[1][2])]; nu=[ppos("neutrinos",k,m) for k,m in enumerate(NU)]
B=[(1,CORN[(1,0,0)],-22,0),(2,np.array([0,0,21.5]),-20,0),(3,(CORN[(1,1,0)]+CORN[(0,1,1)]+CORN[(1,0,1)])/3,60,-40),(4,np.array([0,0,HTOP/2]),24,0),
   (5,np.array([hc*math.tan(52.147*D2R)*math.cos(0.3),hc*math.tan(52.147*D2R)*math.sin(0.3),hc]),22,0),(6,np.array([lr*math.cos(3.6),lr*math.sin(3.6),hc]),-22,0),
   (7,np.array([lr*math.cos(1.57),lr*math.sin(1.57),hc]),0,-20),(8,(lp[1]+lp[2])/2,22,-6),(9,up[2],22,0),(10,nu[1],26,0),(11,(up[1]+dn[1])/2,0,26),
   (12,ax(4.62e22),-24,0),(13,np.zeros(3),-24,6)]
o+=[badge(p,n,dx,dy) for n,p,dx,dy in B]
# twist arc arrow
xa,ya=P(np.array([lr*math.cos(1.57),lr*math.sin(1.57),hc]))
o.append(f'<path d="M{xa-34:.1f},{ya-34:.1f} q34,-16 68,0" fill="none" stroke="#1f2a33" stroke-width="1.8"/><polygon points="{xa+34:.1f},{ya-34:.1f} {xa+22:.1f},{ya-42:.1f} {xa+25:.1f},{ya-29:.1f}" fill="#1f2a33"/>')
y0=70+HG+10
o.append(f'<rect x="20" y="{y0}" width="{W-40}" height="{HL-20}" fill="#ffffff" rx="10" stroke="#d5dde4"/>')
LST=[("1","The cube and its edges: the frame, and the resting lines"),("2","Height: the scale. It orders things, it is not a place in space"),
     ("3","The two pyramids: our side and the far side"),("4","The mirror point at the centre (matter swapped with antimatter)"),
     ("5","Openings: how wide each family's cone is"),("6","The 45° equilibrium (lepton cone, drawn thick)"),
     ("7","The twist: the turn around the axis"),("8","Family ladders: the lines that tie three members together"),
     ("9","Places and link ratios: how far a family sits from equilibrium"),("10","Pure neutrino rungs 1, 2, 3 (only their mixes are ever met)"),
     ("11","The family turn: the tilt between the up and down ladders"),("12","The meeting of the three force lines (traced, never visited)"),
     ("13","The spot and its facing: the perspective")]
for i,(n,t) in enumerate(LST):
    col=i//7; row=i%7; x=44+col*545; yy=y0+38+row*40
    o.append(f'<circle cx="{x+10}" cy="{yy-5}" r="12" fill="#1f2a33"/><text x="{x+10}" y="{yy}" font-size="13" font-weight="700" fill="#fff" text-anchor="middle">{n}</text>')
    o.append(txt(x+32,yy,t,14.8,halo=False))
o.append(txt(44+545+32,y0+38+6*40,"Not drawn, also concepts: the relationship matrix itself, the wave",14.8,"#3b4650",halo=False))
o.append(txt(44,y0+38+7*40+12,"(only landings are met), colour and charge labels, the forces' rules, time's direction, and the 28 measured numbers.",14.8,"#3b4650",halo=False))
o.append(txt(44,y0+38+8*40+22,"A concept is not unreal. It is a relationship. The things in Image 1 are what these relationships produce at their balance points.",15.5,"#1f2a33",700,halo=False))
o.append('</svg>'); open("concepts.svg","w").write("\n".join(o))
