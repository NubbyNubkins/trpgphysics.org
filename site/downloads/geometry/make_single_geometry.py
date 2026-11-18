import numpy as np, math
from scipy.spatial import ConvexHull
n=np.array([1,1,1])/math.sqrt(3)
# frame: z = all-three-together, x = projection of family 1
e1p=np.array([1,0,0])-n/math.sqrt(3); x=e1p/np.linalg.norm(e1p); y=np.cross(n,x)
R=np.vstack([x,y,n])            # world -> (x,y,z) with z vertical
AZ=math.radians(-28); EL=math.radians(16)  # viewing direction
def proj(p):
    q=R@np.asarray(p,float)
    c,s=math.cos(AZ),math.sin(AZ); X=c*q[0]-s*q[1]; Y=s*q[0]+c*q[1]; Z=q[2]
    sx=X; sy=-(Z*math.cos(EL)-Y*math.sin(EL)); depth=Y*math.cos(EL)+Z*math.sin(EL)
    return np.array([sx,sy]),depth
W,H=1400,1520; S=300; OX,OY=500,860
def P(p): v,_=proj(p); return (OX+S*v[0],OY+S*v[1])
def pts(L): return " ".join(f"{a:.1f},{b:.1f}" for a,b in L)
out=[]
add=out.append
add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="Inter, Helvetica, Arial, sans-serif">')
add('<title>TRPG: the single geometry</title>')
add('<defs><linearGradient id="pond" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#eef3f7"/><stop offset="1" stop-color="#dfe8ef"/></linearGradient></defs>')
add(f'<rect width="{W}" height="{H}" fill="url(#pond)"/>')
for i in range(22):   # the pond's fine texture (the grain)
    yy=60+i*68; d="M0,%d "%yy+" ".join(f"Q{x0+35},{yy-5+ (i%2)*2} {x0+70},{yy}" for x0 in range(0,W,70))
    add(f'<path d="{d}" fill="none" stroke="#c9d7e2" stroke-width="1"/>')
# cube on its corner
C=[np.array(v,float) for v in [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]]
edges=[(a,b) for a in range(8) for b in range(a+1,8) if np.sum(np.abs(C[a]-C[b]))==1]
for a,b in edges:
    if a==0: continue
    (x1,y1),(x2,y2)=P(C[a]),P(C[b]); add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#8a97a3" stroke-width="1.6" stroke-dasharray="6 5"/>')
# pyramid (regular tetrahedron in the cube, standing on the bottom corner)
T=[C[0],C[4],C[5],C[6]]
for f in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]:
    add(f'<polygon points="{pts([P(T[i]) for i in f])}" fill="#F5C518" fill-opacity="0.13" stroke="#b08900" stroke-width="1.6" stroke-linejoin="round"/>')
# cones: apex at bottom corner, around the vertical axis
h=0.95
def rim(theta, k=180):
    r=h*math.tan(theta); return [h*n+r*(math.cos(t)*x+math.sin(t)*y) for t in np.linspace(0,2*math.pi,k,endpoint=False)]
cones=[("up quarks",52.147,"#d1495b",3.04),("down quarks",48.150,"#2e8b57",5.71),("leptons",45.0,"#2f6db5",12.733)]
for lab,deg,col,tw in cones:
    R_=rim(math.radians(deg)); P2=[P(p) for p in R_]; allp=np.array([P(C[0])]+P2); hull=ConvexHull(allp)
    add(f'<polygon points="{pts(allp[hull.vertices])}" fill="{col}" fill-opacity="0.07" stroke="none"/>')
    add(f'<polygon points="{pts(P2)}" fill="none" stroke="{col}" stroke-width="2.2"/>')
    for v in hull.vertices:
        if v==0: continue
    # silhouette lines: apex to extreme rim points
    xs=[p[0] for p in P2]; i1,i2=int(np.argmin(xs)),int(np.argmax(xs)); ax,ay=P(C[0])
    for i in (i1,i2): add(f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{P2[i][0]:.1f}" y2="{P2[i][1]:.1f}" stroke="{col}" stroke-width="1.6" stroke-opacity="0.8"/>')
# family edges from the bottom corner (54.74 deg)
fam=[("Family 1",C[1]),("Family 2",C[2]),("Family 3",C[3])]
for lab,v in fam:
    (x1,y1),(x2,y2)=P(C[0]),P(v*1.0); add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#3b4650" stroke-width="3"/>')
    add(f'<circle cx="{x2:.1f}" cy="{y2:.1f}" r="5" fill="#3b4650"/>')
    dx=14 if x2>x1 else -14; anchor="start" if x2>x1 else "end"
    add(f'<text x="{x2+dx:.1f}" y="{y2+5:.1f}" font-size="21" font-weight="600" fill="#3b4650" text-anchor="{anchor}">{lab}</text>')
# vertical axis: all three together
(x1,y1),(x2,y2)=P(C[0]),P(C[7]*1.12)
add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#715800" stroke-width="2" stroke-dasharray="3 6"/>')
add(f'<text x="{x2:.1f}" y="{y2-12:.1f}" font-size="20" fill="#715800" text-anchor="middle" font-weight="600">all three together</text>')
# twist markers on each rim
for lab,deg,col,tw in cones:
    th=math.radians(deg); a=math.radians(tw); r=h*math.tan(th); p=h*n+r*(math.cos(a)*x+math.sin(a)*y)
    (px,py)=P(p); add(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="7" fill="{col}" stroke="#fff" stroke-width="2"/>')
# press direction (19.47 deg)
s=np.array([0.0,1/3,1.0]); s/=np.linalg.norm(s)
d=s[2]*n+s[1]*np.array([1,1,-2])/math.sqrt(6)
(x1,y1),(x2,y2)=P(C[0]),P(d*1.25)
add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#715800" stroke-width="4"/>')
add(f'<circle cx="{x2:.1f}" cy="{y2:.1f}" r="6" fill="#715800"/>')
add(f'<text x="{x2+12:.1f}" y="{y2-6:.1f}" font-size="19" fill="#715800" font-weight="600">the press (19.47°)</text>')
ax,ay=P(C[0]); add(f'<circle cx="{ax:.1f}" cy="{ay:.1f}" r="8" fill="#715800"/>')
add(f'<text x="{ax+16:.1f}" y="{ay+30:.1f}" font-size="19" fill="#715800" font-weight="600">tip: the stable light quark, turning in place by a quarter turn</text>')
# title and legend
add('<text x="40" y="58" font-size="34" font-weight="700" fill="#1f2a33">One geometry</text>')
add('<text x="40" y="92" font-size="19" fill="#3b4650">A cube balanced on one corner, a pyramid standing in it, and nested cones, floating in the still pond</text>')
L=[("#3b4650","Cube edges from the bottom corner: the three families (54.74° from upright)"),
   ("#b08900","Pyramid on its point (fits exactly in the cube): faces lean 19.47°, the press"),
   ("#2f6db5","Leptons: cone at 45.00°, marker turned 12.73° (the twist)"),
   ("#2e8b57","Down quarks: cone at 48.15°, marker turned 5.71°"),
   ("#d1495b","Up quarks: cone at 52.15°, marker turned 3.04°"),
   ("#715800","Upright line: all three together (equal weight in every family)")]
y0=140
for col,t in L:
    add(f'<rect x="1000" y="{y0-14}" width="22" height="14" rx="3" fill="{col}"/>'); 
    words=t.split(); line="";lines=[]
    for w in words:
        if len(line+w)>34: lines.append(line); line=""
        line+=w+" "
    lines.append(line)
    for j,ln in enumerate(lines): add(f'<text x="1032" y="{y0+j*22:.0f}" font-size="17" fill="#1f2a33">{ln.strip()}</text>')
    y0+=22*len(lines)+16
add(f'<text x="1000" y="{y0+10}" font-size="16" fill="#3b4650">The closer a cone sits to the cube</text>')
add(f'<text x="1000" y="{y0+31}" font-size="16" fill="#3b4650">edges, the more lopsided that</text>')
add(f'<text x="1000" y="{y0+52}" font-size="16" fill="#3b4650">family\'s weights. The two quark</text>')
add(f'<text x="1000" y="{y0+73}" font-size="16" fill="#3b4650">ladders not sharing one frame is</text>')
add(f'<text x="1000" y="{y0+94}" font-size="16" fill="#3b4650">the quark mixing triangle.</text>')
add(f'<text x="1000" y="{y0+130}" font-size="16" fill="#3b4650">Background ripples: the pond and</text>')
add(f'<text x="1000" y="{y0+151}" font-size="16" fill="#3b4650">its fine grain. The far side,</text>')
add(f'<text x="1000" y="{y0+172}" font-size="16" fill="#3b4650">neutrino waves and gravity act</text>')
add(f'<text x="1000" y="{y0+193}" font-size="16" fill="#3b4650">from the water, not inside the cube.</text>')
for j,ln in enumerate(['The three shapes found in every dimension:','triangle family = the pyramid (tetrahedron);','square family = the cube;','cross family = the octahedron, whose six','corners are the families and their opposites.']):
    add(f'<text x="1000" y="{y0+230+j*21}" font-size="16" fill="#715800" font-weight="{600 if j==0 else 400}">{ln}</text>')

# ---------- lower panel: the two quark ladders as opposite sides of one cone ----------
add('<line x1="40" y1="1010" x2="1360" y2="1010" stroke="#b9c7d3" stroke-width="1.5"/>')
add('<text x="40" y="1050" font-size="26" font-weight="700" fill="#1f2a33">The two quark ladders: opposite sides of one cone</text>')
add('<text x="40" y="1078" font-size="17" fill="#3b4650">Up-type and down-type are the one quark pair seen two ways. Height = family weight; width = how many times heavier the up-type partner is.</text>')
mu_=[1.242,627.8,162500.0]; md_=[2.678,53.36,2873.0]; nu_=["up","charm","top"]; nd_=["down","strange","bottom"]
CX=330; KX=88; YT_,YB_=1110,1440
def yv(m): return YB_-(math.log10(m)-0.0)/(4.6)*(YB_-YT_)
U_=[];D_=[]
for k in range(3):
    gm=math.sqrt(mu_[k]*md_[k]); w=math.log(mu_[k]/md_[k]); yy=yv(gm)
    U_.append((CX+KX*w/2,yy)); D_.append((CX-KX*w/2,yy))
# centre line (all three together)
add(f'<line x1="{CX}" y1="{YT_-10}" x2="{CX}" y2="{YB_+10}" stroke="#715800" stroke-width="1.5" stroke-dasharray="3 6"/>')
# crossing point
f=0.238; yc=U_[0][1]+f*(U_[1][1]-U_[0][1])
add(f'<polyline points="{pts(U_)}" fill="none" stroke="#d1495b" stroke-width="3"/>')
add(f'<polyline points="{pts(D_)}" fill="none" stroke="#2e8b57" stroke-width="3"/>')
add(f'<circle cx="{CX}" cy="{yc:.1f}" r="7" fill="#715800"/>')
add(f'<text x="{CX-14}" y="{yc+5:.1f}" font-size="16" fill="#715800" font-weight="600" text-anchor="end">sides cross (the tip)</text>')
for k in range(3):
    (ux,uy),(dx_,dy)=U_[k],D_[k]
    add(f'<line x1="{dx_:.1f}" y1="{dy:.1f}" x2="{ux:.1f}" y2="{uy:.1f}" stroke="#8a97a3" stroke-width="1.2" stroke-dasharray="4 4"/>')
    add(f'<circle cx="{ux:.1f}" cy="{uy:.1f}" r="7" fill="#d1495b" stroke="#fff" stroke-width="2"/>')
    add(f'<circle cx="{dx_:.1f}" cy="{dy:.1f}" r="7" fill="#2e8b57" stroke="#fff" stroke-width="2"/>')
    au="start" if ux>=dx_ else "end"; ad="end" if ux>=dx_ else "start"; ou=12 if ux>=dx_ else -12
    add(f'<text x="{ux+ou:.1f}" y="{uy+5:.1f}" font-size="16" fill="#d1495b" font-weight="600" text-anchor="{au}">{nu_[k]}</text>')
    add(f'<text x="{dx_-ou:.1f}" y="{dy+5:.1f}" font-size="16" fill="#2e8b57" font-weight="600" text-anchor="{ad}">{nd_[k]}</text>')
rat=["×0.46, crossed","×12","×57"]
for k in range(3):
    add(f'<text x="600" y="{U_[k][1]+5:.1f}" font-size="16" fill="#3b4650">family {k+1}: up-type {rat[k]}</text>')
# right-hand notes
NX=865; ny=1120
notes=[("#1f2a33",600,"What the picture shows"),
 ("#3b4650",400,"The sides meet just above the first family,"),
 ("#3b4650",400,"where up and down swap places (down is the"),
 ("#3b4650",400,"heavier one there). That is the tip above."),
 ("#3b4650",400,"The V opens upward, but slower each step:"),
 ("#3b4650",400,"×25, then ×5. Curved like a tulip, widest at top."),
 ("#2e8b57",400,"Down rungs: close at the bottom (×20), wider at the top (×54)."),
 ("#d1495b",400,"Up rungs: wide at the bottom (×506), closer at the top (×259)."),
 ("#2f6db5",400,"Leptons, for comparison: ×207, then ×17."),
 ("#1f2a33",600,"How fast it opens (second step ÷ first, as powers)"),
 ("#3b4650",400,"0.49 at our energies, rising to 0.53 (the lepton"),
 ("#3b4650",400,"ladder's own rate) near 2×10⁸ GeV, where the"),
 ("#3b4650",400,"far-side weak pairs sit. Any rate from 0.49 to 0.55"),
 ("#3b4650",400,"matches at some energy, so this is a hint, not a rule.")]
for col,wt,t in notes:
    if wt==600 and ny>1120: ny+=10
    add(f'<text x="{NX}" y="{ny}" font-size="16" fill="{col}" font-weight="{wt}">{t}</text>'); ny+=23

add('<text x="40" y="1500" font-size="14" fill="#6b7782">© 2026 Matthew Gelfman · trpgphysics.org · angles computed exactly (lepton cone: Koide 1982, Brannen 2006)</text>')
add('</svg>')
open('geometry/trpg_single_geometry.svg','w').write("\n".join(out))
print("ok")
