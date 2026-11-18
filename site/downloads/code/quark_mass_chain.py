"""Quark-mass chain (14 Nov 2026): one input (bottom mass) + three relations already in the matrices.
1) top fully linked to the Higgs: y_t = 1 -> m_t = v/sqrt2  (near-fixed-point: Pendleton & Ross 1981; Hill 1981)
2) heavy triples obey the 45-degree balance (Koide Q = 2/3): (c,b,t) then (s,c,b) with the strange root negative (Rivero 2011)
3) Cabibbo from masses: m_d = m_s |Vus|^2 (Gatto, Sartori & Tonin 1968)
Caveat: masses mix schemes (top pole; c, b at own scale; s, d at 2 GeV)."""
import math
def koide_light(m1,m2,sign=1):
    a,b=math.sqrt(m1),math.sqrt(m2); B=-4*sign*(a+b); C=3*(m1+m2)-2*(a+b)**2
    return [x*x for x in [(-B-math.sqrt(B*B-4*C))/2,(-B+math.sqrt(B*B-4*C))/2] if x>0]
mb=4183.0; mt=246219.65/math.sqrt(2)
mc=min(koide_light(mt,mb)); ms=min(x for x in koide_light(mb,mc,-1) if x<mc); md=ms*0.22431**2
for n,p,m in [("top",mt,172570),("charm",mc,1273),("strange",ms,93.5),("down",md,4.70)]:
    print(f"{n:8s} predicted {p:10.2f}  measured {m:10.2f}  off {100*(p-m)/m:+.1f}%")
print("up: no relation reaches it (open)")

# --- an earlier test additions: the up quark through the meson cells (Weinberg 1977), and the blind spot of the lean ---
Kp,K0,pp,p0=493.677,497.611,139.57039,134.9768
r=(Kp**2-K0**2+2*p0**2-pp**2)/(K0**2-Kp**2+pp**2)   # m_u/m_d from kaon and pion cells, leading order
print(f"up/down from meson cells (Weinberg 1977): {r:.3f}; lattice with electromagnetic corrections: 0.465")
print(f"up predicted = m_d x ratio: {md*r:.2f} MeV (leading order), {md*0.465:.2f} MeV (lattice ratio); measured 2.16")
def ang(ms):
    q=sum(ms)/sum(math.sqrt(m) for m in ms)**2; return math.degrees(math.acos(math.sqrt(1/(3*q))))
print("up-family angle with m_u = 0, 2.16, 10 MeV:", [round(ang((u,1273,172570)),3) for u in (0,2.16,10)], "(blind to the up mass)")

# --- an earlier test: up quark WITHOUT the lattice electric correction ---
# eta -> 3 pions measures the up-down mass difference with electric effects suppressed (Sutherland 1966; Bell & Sutherland 1968).
# Leutwyler's ellipse: (m_u/m_d)^2 + (m_s/m_d)^2 / Q^2 = 1 ; Q = 22.1 +/- 0.7 (Colangelo, Lanz, Leutwyler & Passemar 2018, from eta->3pi data)
rsd=1/0.22431**2   # m_s/m_d from the Cabibbo cell (our chain)
for Qv in (21.4,22.1,22.8):
    x=math.sqrt(1-(rsd/Qv)**2); print(f"Q={Qv}: m_u/m_d={x:.3f}, up={md*x:.2f} MeV")
