from onepic import *
F_='font-family="Inter, Helvetica, Arial, sans-serif"'
hM=hgt(2.41e17); hF=hgt(4.56e22); hG=hgt(8.63e27); hD=hF
ex=[]; lines=[]
for k,ph in enumerate([30,150,270]):
    a=ph*D2R; rr=0.38*rcube(hM,a); p=np.array([rr*math.cos(a),rr*math.sin(a),hM])
    ex.append(dict(p=p,r=11,fill="#e67e22",stroke="#7a3e00",label=f"MMG{k+1}" if k<3 else ""))
    rd=0.30*rcube(hD,a); pd=np.array([rd*math.cos(a),rd*math.sin(a),hD])
    ex.append(dict(p=pd,r=8,fill="#c0392b",stroke="#6e1e16",label=f"D{k+1}" if k==0 else ""))
    lines.append((p,pd,"#e67e22","4 3"))
ax=lambda h: np.array([0,0,h])
lines.append((ax(hG),ax(hF),"#8e3b8e","8 5")); lines.append((ax(hF),ax(hM),"#8e3b8e","8 5"))
ex.append(dict(p=ax(hF),r=7,fill="#fff",stroke="#c0392b",label="",lc="#c0392b"))
W,H=1300,1480
s=render((-35,22),W,1300,mirror=True,marks=True,labels=True,extra=ex,extra_lines=lines)
s=s.replace(f'viewBox="0 0 {W} 1300" width="{W}" height="1300"',f'viewBox="0 0 {W} {H}" width="{W}" height="{H}"').replace(f'<rect width="{W}" height="1300" fill="#eef3f7"/>',f'<rect width="{W}" height="{H}" fill="#eef3f7"/>')
head=(f'<rect x="0" y="0" width="{W}" height="92" fill="#eef3f7" fill-opacity="0.94"/>'
 f'<text x="28" y="40" font-size="28" font-weight="700" fill="#1f2a33" {F_}>Where the MMG sits in the one picture</text>'
 f'<text x="28" y="72" font-size="17" fill="#3b4650" {F_}>Three MMGs (orange), one per family, on the turned-over (hidden) half; their coloured partners D (red) higher up at the meeting point.</text>')
rows=["Orange: the three MMGs, 2.4×10⁸ GeV each. Weak charge, no colour, spin ½, charges 0 and ±1. Both-handed (vector-like).",
      "Red: the coloured partners D at the meeting point, 4.6×10¹³ GeV (2×10⁵ heavier). Orange dashes: each MMG with its own D.",
      "Purple dashes on the centre line: the mirror. The grain sits as far above the meeting point as the MMGs sit below it.",
      "They point the opposite way to our three families (30°, 150°, 270° against 210°, 330°, 90°): the turned-over half.",
      "Where the strengths change: above the MMG height the weak line bends, so all three lines meet at the meeting point."]
L=[f'<text x="40" y="{1250+i*34}" font-size="16" fill="#1f2a33" {F_}>{t}</text>' for i,t in enumerate(rows)]
foot=f'<text x="28" y="{H-18}" font-size="13" fill="#5b6670" {F_}>© 2026 Matthew Gelfman · trpgphysics.org · weights from the force strengths (two loop) and the grain mirror</text>'
s=s.replace('</svg>',head+"".join(L)+foot+'</svg>')
open("mmg.svg","w").write(s)
