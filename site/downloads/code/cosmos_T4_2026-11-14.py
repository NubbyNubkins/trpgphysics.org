"""TRPG and the cosmos (test T4), 14 Nov 2026. Planck 2018: H0=67.4, Omega_L=0.685, Omega_k=0.001+-0.002."""
import math
c=2.99792458e8; G=6.674e-11; hbarc_MeVfm=197.3269804; Mpc=3.0857e22; MeV_J=1.602176634e-13
H0=67.4e3/Mpc; rho_c=3*H0**2/(8*math.pi*G)              # kg/m^3
uL=0.685*rho_c*c**2                                      # J/m^3 dark energy density
uL_MeVfm3=uL/MeV_J*1e-45
B=66.0
print("dark energy density = %.3e MeV/fm^3 ; proton-scale baseline B = %.0f MeV/fm^3 ; ratio B/rho_L = %.1e"%(uL_MeVfm3,B,B/uL_MeVfm3))
# energy scales (fourth roots) in natural units: E = (u * (hbar c)^3)^(1/4)
EL=(uL_MeVfm3*hbarc_MeVfm**3)**0.25*1e9   # meV
EB=(B*hbarc_MeVfm**3)**0.25
print("baseline energy scale: proton layer %.0f MeV ; dark-energy layer %.2f meV"%(EB,EL))
print("compare M2 neutrino masses (meV): 0.36, 8.7, 50.1 ; sqrt(m1*m2) = %.2f meV ; m2/4 = %.2f"%(math.sqrt(0.36*8.66),8.66/4))
# Hubble-sphere compactness
R=c/H0; M=(4*math.pi/3)*R**3*rho_c
print("Hubble sphere: R = %.2e m, M = %.2e kg, compactness 2GM/(Rc^2) = %.4f (exactly 1 for a flat universe)"%(R,M,2*G*M/(R*c**2)))
