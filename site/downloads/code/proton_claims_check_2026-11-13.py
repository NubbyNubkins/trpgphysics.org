"""Re-run of the earlier version's proton 'predictions' using its own formulas and typed constants (13 Nov 2026)."""
import numpy as np
MN,MD,mpi=1085-3*48.91,1085+3*48.91,139.57      # earlier version: centroid 1085 & C_ss 48.91 typed -> split = 6*48.91 = measured gap
p=np.sqrt((MD**2-(MN+mpi)**2)*(MD**2-(MN-mpi)**2))/(2*MD); EN=np.sqrt(MN**2+p**2)
G=lambda f:(f**2/(12*np.pi*mpi**2))*((EN+MN)/(2*MD))*p**3
print("Delta width with the earlier version f=2.18: %.1f MeV (the earlier version claimed 114.2; measured 117). f needed for 117: %.2f"%(G(2.18),2.18*np.sqrt(117/G(2.18))))
Z2=1/1.214
print("gA from stated formula: %.4f (the earlier version then adds +0.0517 by hand -> 1.2781)"%(Z2*(1.5129-(2/27)*0.214*1.5129)))
print("mu_p, mu_n from Z2*bare: %.3f, %.3f (the earlier version typed 2.368, -1.518, then +/-0.412)"%(Z2*2.5424,Z2*-1.6949))
print("Delta m_np: uses PDG m_d-m_u as input and multiplies EM energy by a chosen 2.85")
print("von Laue: wall tension rescaled by scale_b=-int_q/int_b, so residual = 0 by construction")
