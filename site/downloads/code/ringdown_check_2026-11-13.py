"""Black-hole ringdown and echo test for TRPG's 'central particle' black hole (13 Nov 2026).
Kerr l=m=2 fundamental mode fits: Berti, Cardoso & Will (2006). Echo delay: Abedi, Dykaar & Afshordi (2017).
GW150914 remnant: source mass 62 Msun, spin 0.67, z = 0.09 (LIGO/Virgo 2016)."""
import math
Msun_s=4.925491e-6; Msun_m=1476.625
M=62*(1+0.09); a=0.67
Mw=1.5251-1.1568*(1-a)**0.1292; Q=0.7000+1.4187*(1-a)**(-0.4990)
w=Mw/(M*Msun_s); f=w/(2*math.pi); tau=2*Q/w
print("GR (Kerr) ringdown: f = %.0f Hz, tau = %.1f ms   | measured 251 +/- 8 Hz, 4.0 +/- 0.3 ms"%(f,tau*1e3))
lP=1.616255e-35
dt=4*M*Msun_s*(1+1/math.sqrt(1-a*a))*math.log(M*Msun_m/lP)
print("If a reflecting central-particle surface sat a Planck length above the horizon: echoes every %.3f s"%dt)
# Compactness ladder C = 2GM/(Rc^2)
for n,C in [("Sun",4.24e-6),("white dwarf (Sirius B)",2*1.02*Msun_m/5.8e6),("neutron star 1.4 Msun, R 11.4 km (K1)",2*1.4*Msun_m/11.4e3),
            ("Buchdahl limit, any static star",8/9),("black hole horizon",1.0)]:
    print("%-40s C = %.3g"%(n,C))
