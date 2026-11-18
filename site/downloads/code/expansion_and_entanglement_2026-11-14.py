"""M11 expansion from the 'full bag / full glass' and M12 uncertainty as wave math. 14 Nov 2026.
Planck 2018 values. Friedmann (1922), Lemaitre (1927)."""
import numpy as np, math
from scipy.integrate import quad
H0=67.4; Om=0.315; Or=9.1e-5; OL=1-Om-Or   # full glass: parts sum to exactly 1
print("full glass: matter %.3f + dark energy %.3f + radiation %.5f = %.4f (Planck: curvature 0.001 +/- 0.002)"%(Om,OL,Or,Om+OL+Or))
H=lambda a: H0*math.sqrt(Om/a**3+Or/a**4+OL)
Mpc_km=3.0857e19; Gyr_s=3.156e16
age=quad(lambda a: 1/(a*H(a)),1e-10,1)[0]*Mpc_km/Gyr_s
print("age from the expansion equation: %.2f billion years (Planck 2018: 13.787 +/- 0.020)"%age)
for a in [1e-9,1e-4,0.5,1,2]:
    print("size %.0e of today: matter share %.3f, radiation %.3f, dark energy %.3f"%(a,(Om/a**3)/(H(a)/H0)**2,(Or/a**4)/(H(a)/H0)**2,OL/(H(a)/H0)**2))
# local conservation check: rho_m * a^3 constant (same amount of 'water'), rho_L constant per volume
print("matter amount in a stretching region (rho*a^3): constant -> same water; dark energy per volume constant -> its total grows: global total not conserved, only local (as in Einstein's theory)")
# M12: uncertainty as pure wave math (Fourier): Gaussian pulse gives dx*dk = 1/2 exactly
x=np.linspace(-50,50,2**16); s=2.0; psi=np.exp(-x**2/(4*s*s)); p=abs(psi)**2; p/=np.trapezoid(p,x)
dx=math.sqrt(np.trapezoid(x**2*p,x))
k=np.fft.fftshift(np.fft.fftfreq(len(x),x[1]-x[0]))*2*np.pi; P=abs(np.fft.fftshift(np.fft.fft(psi)))**2; P/=np.trapezoid(P,k)
dk=math.sqrt(np.trapezoid(k**2*P,k)); print("wave pulse: spread in place x spread in wavenumber = %.4f (Heisenberg minimum 0.5)"%(dx*dk))
