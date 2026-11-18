# When and where on the cone each kind of particle drops out of the stretching light
import numpy as np
from scipy.integrate import quad
lg=np.log10
GRAIN,H0eV=8.633e27,1.438e-33
T0=2.349e-4            # eV, light today (2.725 K)
Tnu0=1.676e-4          # eV, relic neutrino temperature today (1.945 K)
H0=67.4/3.0857e19; Om=0.315; Orad=9.15e-5; OL=1-Om-Orad; yr=3.156e7
def gstar(T):          # T in eV, standard counts
    for lim,g in ((2e11,106.75),(8e10,96),(4e9,86),(1.2e9,75),(2e8,62),(1.5e8,17.25),(1e8,14),(5e5,10.75)):
        if T>lim: return g
    return 3.36
def age(T):            # seconds since the start when the light's energy (kT) was T
    if T>10:           # radiation era
        return 2.42/np.sqrt(gstar(T))/(T/1e6)**2
    z=T/T0-1; f=lambda a: 1/(a*H0*np.sqrt(Orad/a**4+Om/a**3+OL)); return quad(f,0,1/(1+z))[0]
def fmt(t):
    if t<1e-3: return "%.1e s"%t
    if t<120: return "%.2g s"%t
    if t<3*3600: return "%.0f min"%(t/60)
    if t<yr: return "%.2g days"%(t/86400)
    y=t/yr; return ("%.3g thousand years"%(y/1e3) if y<1e6 else ("%.3g million years"%(y/1e6) if y<1e9 else "%.3g billion years"%(y/1e9)))
nu=[3.6e-4,8.62e-3,5.01e-2]; le=[0.511e6,105.66e6,1776.9e6]; N=[le[i]**2/nu[i] for i in range(3)]
ev=[("the grain: start",GRAIN),("D, coloured far-side parts",4.62e22),("hidden partner 3",N[2]),("hidden partner 2",N[1]),("MMG",2.41e17),("hidden partner 1",N[0]),
    ("top",172.6e9),("the surface's resistance switches on (Higgs settles)",160e9),("Higgs",125.2e9),("Z",91.2e9),("W",80.4e9),("bottom",4.18e9),("tau",1.777e9),("charm",1.27e9),
    ("quarks lock into protons and neutrons",155e6),("muon",105.7e6),("strange",93.5e6),("neutrinos stop touching anything",1.0e6),("electron (pairs go, the leftover stays)",0.511e6),
    ("first nuclei",7e4),("majoron (dark matter candidate)",8.8e3),("atoms form; light set free",0.26),("neutrino 3 slows (stops behaving like light)",None),("neutrino 2 slows",None),("neutrino 1 slows",None),("TODAY",T0)]
print("(1) Stretch check: grain-length wave carried by expansion alone to today")
s=(GRAIN/T0)*(106.75/3.91)**(1/3); print("   stretch factor %.2e ; 1.436e-34 m -> %.3g m ; measured glow: typical wave %.3g m"%(s,1.436e-34*s,1.2398e-6/(2.7*T0)))
print("\n(2) When the falling light passed each level (light's energy = that particle's energy)")
for n,E in ev:
    if E is None:
        m={"neutrino 3 slows (stops behaving like light)":nu[2],"neutrino 2 slows":nu[1],"neutrino 1 slows":nu[0]}[n]; z=m/(3.15*Tnu0)-1; T=T0*(1+z)
        print(f"   {n:55s} light at {T:.3g} eV   stretch since then x{1+z:.3g}   age {fmt(age(T))}")
    else: print(f"   {n:55s} light at {E:.3g} eV   stretch since then x{(E/T0)*(gstar(E)/3.91)**(1/3):.3g}   age {fmt(age(E))}")
print("\n(3) How much of each relic thing is free movement (stretches) and how much is held (does not), today")
p=3.15*Tnu0
for n,m in (("light",0),("neutrino 1",nu[0]),("neutrino 2",nu[1]),("neutrino 3",nu[2])):
    E=np.sqrt(p*p+m*m) if m else 1; print(f"   {n:12s} moving part {100*((E-m)/E if m else 1):.4g}%   held part {100*(m/E if m else 0):.4g}%")
Te=T0  # a free electron at today's light temperature
print("   electron     moving part %.2g%%"%(100*1.5*Te/0.511e6)); print("   proton       moving part %.2g%%"%(100*1.5*Te/938.27e6))
print("\n(4) The floor of the stretch: the surface's own glow (Gibbons & Hawking 1977) = surface energy x sqrt(dark share)/2pi = %.3g eV"%(H0eV*np.sqrt(OL)/(2*np.pi)))
print("    light today is at %.3g eV, %.1f tens above the surface; the cone is %.1f tens tall; light has come %.1f%% of the way down"%(T0,lg(T0/H0eV),lg(GRAIN/H0eV),100*lg(GRAIN/T0)/lg(GRAIN/H0eV)))
