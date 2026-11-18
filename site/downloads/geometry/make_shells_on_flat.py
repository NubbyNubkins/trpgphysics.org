import math, numpy as np
F_='font-family="Inter, Helvetica, Arial, sans-serif"'
W,H=1150,1530; CX,CY=470,600; K=1.0
def hg(eV):
    h=math.log10(eV)+5
    return 45+(h-1.5)*14 if h<8 else 136+(h-8)*37
o=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" {F_}><rect width="{W}" height="{H}" fill="#f4f7fa"/>']
def txt(x,y,s,sz=15,c="#1f2a33",w=400,a="start",halo=True):
    h=' paint-order="stroke" stroke="#f4f7fa" stroke-width="4"' if halo else ''
    o.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{sz}" fill="{c}" font-weight="{w}" text-anchor="{a}"{h}>{s}</text>')
txt(30,44,"Where the shells appear on the flat image",27,w=700,halo=False)
txt(30,72,"Distance from the centre shows scale: further out means more energy and a smaller size. The three directions are the three rungs.",15,"#3b4650",halo=False)
RM=hg(172.57e9)*K+22
o.append(f'<circle cx="{CX}" cy="{CY}" r="{RM:.1f}" fill="#e6eef5" stroke="#9fb3c4" stroke-width="2"/>')
def pol(r,deg): a=math.radians(deg); return CX+r*math.cos(a),CY-r*math.sin(a)
DEG=[210,330,90]
for d in DEG:
    x,y=pol(RM-4,d); o.append(f'<line x1="{CX}" y1="{CY}" x2="{x:.1f}" y2="{y:.1f}" stroke="#9fb3c4" stroke-width="1.2" stroke-dasharray="3 5"/>')
# shells
GOLD="#b08900"; DG="#715800"
ra=hg(3.73e3)*K; rn1=hg(27e6)*K; rn2=hg(164e6)*K; rp=hg(235e6)*K
o.append(f'<circle cx="{CX}" cy="{CY}" r="{(rn1+rn2)/2:.1f}" fill="none" stroke="{GOLD}" stroke-opacity="0.28" stroke-width="{rn2-rn1:.1f}"/>')
o.append(f'<circle cx="{CX}" cy="{CY}" r="{ra:.1f}" fill="none" stroke="{DG}" stroke-width="3.5"/>')
o.append(f'<circle cx="{CX}" cy="{CY}" r="{rp:.1f}" fill="none" stroke="{DG}" stroke-width="3.5"/>')
NU=[3.6e-4,math.sqrt(3.6e-4**2+7.42e-5),math.sqrt(3.6e-4**2+2.51e-3)]
FAM=[("neutrinos","#7a5aa6",[("",NU[0]),("",NU[1]),("",NU[2])]),
     ("leptons","#2f6db5",[("electron",0.511e6),("muon",105.66e6),("tau",1776.9e6)]),
     ("down","#2e8b57",[("down",4.70e6),("strange",93.5e6),("bottom",4183e6)]),
     ("up","#d1495b",[("up",2.16e6),("charm",1273e6),("top",172570e6)])]
LASTS={"electron","up","down"}
OFF={"electron":(4,-16,"middle"),"up":(-16,6,"end"),"down":(10,24,"start"),"muon":(8,-16,"middle"),"strange":(-14,22,"end"),"charm":(14,4,"start"),
     "tau":(14,10,"start"),"bottom":(-14,4,"end"),"top":(16,5,"start")}
for name,col,it in FAM:
    P=[pol(hg(m)*K,d) for (_,m),d in zip(it,DEG)]
    o.append(f'<polygon points="{" ".join(f"{a:.1f},{b:.1f}" for a,b in P)}" fill="{col}" fill-opacity="0.07" stroke="{col}" stroke-width="1.8" stroke-opacity="0.8"/>')
for name,col,it in FAM:
    for (lab,m),d in zip(it,DEG):
        x,y=pol(hg(m)*K,d)
        if lab in LASTS: o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="9" fill="{col}" stroke="#1f2a33" stroke-width="2.2"/>')
        else: o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="#fff" stroke="{col}" stroke-width="2"/>')
        if lab: dx,dy,a=OFF[lab]; txt(x+dx,y+dy,lab,14,col,600,a)
o.append(f'<circle cx="{CX}" cy="{CY}" r="10" fill="#e0a800" stroke="#1f2a33" stroke-width="2"/>'); txt(CX+14,CY-12,"light",14,"#8a6d00",600)
txt(CX+62,CY+62,"neutrinos",13.5,"#7a5aa6",600)
# W Z H
x,y=pol(hg(1.0e11)*K,150); 
o.append(f'<circle cx="{CX}" cy="{CY}" r="{hg(1.0e11)*K:.1f}" fill="none" stroke="#715800" stroke-width="1.2" stroke-dasharray="2 6" stroke-opacity="0.7"/>'); txt(x-6,y-8,"W, Z, Higgs",13.5,"#715800",600,"end")
# ring labels (to the right, along 20 deg)
for r,t1,t2 in ((ra,"ATOM shell","held by the electric pull (light)"),((rn1+rn2)/2,"NUCLEUS shell","held by the leftover strong pull"),(rp,"PROTON and NEUTRON shell","held by the strong pull (gluons)")):
    pass
def lab(r,deg,t1,t2,ax,ay,anc="start"):
    x,y=pol(r,deg); o.append(f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{ax}" y2="{ay}" stroke="{DG}" stroke-width="1.3"/>')
    txt(ax+(6 if anc=="start" else -6),ay-2,t1,15,DG,700,anc); txt(ax+(6 if anc=="start" else -6),ay+16,t2,13.5,"#3b4650",400,anc)
lab(ra,35,"ATOM shell","held by the electric pull",935,330)
lab((rn1+rn2)/2,8,"NUCLEUS shell (band)","held by the leftover",935,500)
txt(941,534,"strong pull",13.5,"#3b4650")
lab(rp,-14,"PROTON and NEUTRON","shell (one ring, same size)",935,700)
txt(941,734,"held by the strong pull",13.5,"#3b4650")
txt(CX,CY-RM-12,"third rung",14,"#5b6b78",700,"middle"); txt(60,CY+RM*0.55+60,"first rung (what lasts)",15,"#1f2a33",700); txt(CX+RM*0.80,CY+RM*0.55+60,"second rung",14,"#5b6b78",700,"middle")
txt(CX-70,CY-88,"empty stretch",12.5,"#8a97a3",400,"end"); txt(CX-70,CY-73,"(drawn shortened)",12.5,"#8a97a3",400,"end")
y0=1105
o.append(f'<rect x="20" y="{y0}" width="{W-40}" height="{H-y0-20}" fill="#ffffff" rx="10" stroke="#d5dde4"/>')
txt(40,y0+34,"What the picture shows",20,w=700,halo=False)
L=[("b","Three rings, in order outward: atom, nucleus, then proton and neutron (those two share one ring)."),
   ("","Each ring is the balance point of one pull: electric for the atom, leftover strong for the nucleus, strong for the proton and neutron."),
   ("","The weak pull and the Higgs link make no ring at all. They bind nothing."),
   ("b","The rings close around the first rung only. Electron, up and down are the only filled dots, and they are the only ones inside a lasting shell."),
   ("","Strange and the muon sit right at the nucleus ring. They can get inside a shell for a moment (strange particles, muon atoms) but do not last."),
   ("","Charm, tau, bottom, W, Z, Higgs and top lie outside every ring. The top breaks down before a shell can even form."),
   ("b","The triangles do widen outward, but they are lopsided and cross each other. They do not each get their own ring."),
   ("","In real space a shell wraps each single bubble. On this chart the rings mark the scale at which that wrapping happens.")]
for i,(b,s) in enumerate(L): txt(44,y0+68+i*40,s,14.8,"#1f2a33",700 if b else 400,halo=False)
o.append('</svg>'); open("rings.svg","w").write("\n".join(o))
print("atom",3.73e3,"nucleon",197.327/0.84,"MeV; top life 5e-25 s vs shell time ~3e-24 s")
