from onepic import *
W,H=1300,1480
s=render((-35,22),W,1300,mirror=True,marks=True,decay='fwd',labels=True)
s=s.replace(f'viewBox="0 0 {W} 1300" width="{W}" height="1300"',f'viewBox="0 0 {W} {H}" width="{W}" height="{H}"')
s=s.replace(f'<rect width="{W}" height="1300" fill="#eef3f7"/>',f'<rect width="{W}" height="{H}" fill="#eef3f7"/>')
F='font-family="Inter, Helvetica, Arial, sans-serif"'
head=(f'<rect x="0" y="0" width="{W}" height="92" fill="#eef3f7" fill-opacity="0.94"/>'
 f'<text x="28" y="40" font-size="28" font-weight="700" fill="#1f2a33" {F}>The one picture, with everything in it</text>'
 f'<text x="28" y="72" font-size="17" fill="#3b4650" {F}>Weights, waves, mixing, spin, charge, colour, the mirror half, and the way things decay as time runs forward.</text>')
leg=[("#c9a227","solid","gold pyramid: our side of the bag (the matrix seen from one spot)"),
     ("#2b7bb9","dash","blue pyramid: the mirror half (antimatter, the far side turned around)"),
     ("#16a085","dash","teal lines: the octahedron where the two pyramids overlap (Kepler's star, 1619)"),
     ("#c0392b","arrow","red arrows: decays as time runs forward (heavy settles toward light)"),
     ("#ffffff","circ","white circles: the blends the weak force meets (neutrino and quark mixing)"),
     ("#3b4650","text","marks in each bubble: half-arc = spin ½, number = charge, three dots = three colours (quarks only)")]
y0=1250; L=[]
for i,(c,k,t) in enumerate(leg):
    x=40; y=y0+i*34
    if k=="solid": L.append(f'<line x1="{x}" y1="{y}" x2="{x+40}" y2="{y}" stroke="{c}" stroke-width="4"/>')
    elif k=="dash": L.append(f'<line x1="{x}" y1="{y}" x2="{x+40}" y2="{y}" stroke="{c}" stroke-width="3" stroke-dasharray="7 4"/>')
    elif k=="arrow": L.append(f'<line x1="{x}" y1="{y}" x2="{x+32}" y2="{y}" stroke="{c}" stroke-width="2.2"/><polygon points="{x+40},{y} {x+31},{y-5} {x+31},{y+5}" fill="{c}"/>')
    elif k=="circ": L.append(f'<circle cx="{x+20}" cy="{y}" r="9" fill="#fff" stroke="#2e8b57" stroke-width="2"/>')
    else: L.append(f'<text x="{x+20}" y="{y+6}" font-size="18" text-anchor="middle" fill="{c}" {F}>½ ⅔</text>')
    L.append(f'<text x="{x+54}" y="{y+6}" font-size="16" fill="#1f2a33" {F}>{t}</text>')
foot=f'<text x="28" y="{H-18}" font-size="13" fill="#5b6670" {F}>© 2026 Matthew Gelfman · trpgphysics.org · masses: Particle Data Group 2024 · mixing: NuFIT 2024</text>'
s=s.replace('</svg>',head+"".join(L)+foot+'</svg>')
open("full_leg.svg","w").write(s)
