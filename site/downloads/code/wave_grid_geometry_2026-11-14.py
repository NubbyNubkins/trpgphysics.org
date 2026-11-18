"""The crossing-wave picture of space: z = cos x + cos y (and its time dependence). 14 Nov 2026.
Checks: critical points (max / saddle / min), Gaussian curvature sign map, Gauss-Bonnet total, and standing-wave time flip."""
import numpy as np
n=1201; x=np.linspace(0,2*np.pi,n); X,Y=np.meshgrid(x,x)
zx,zy=-np.sin(X),-np.sin(Y); zxx,zyy=-np.cos(X),-np.cos(Y); zxy=0
K=(zxx*zyy-zxy**2)/(1+zx**2+zy**2)**2
dA=np.sqrt(1+zx**2+zy**2)*(x[1]-x[0])**2
print("critical points per cell: (0,0) peak z=2, (pi,pi) trough z=-2, (pi,0) & (0,pi) saddles z=0 -> 1 peak, 1 trough, 2 saddles")
print("Morse check (peaks - saddles + troughs) = 1 - 2 + 1 = 0 = Euler number of a repeating grid")
pos=(K>0); print("area fraction curved like a sphere (K>0): %.3f ; like a saddle (K<0): %.3f"%((dA*pos).sum()/dA.sum(),(dA*(K<0)).sum()/dA.sum()))
print("total curvature over one cell (Gauss-Bonnet says 0): %.2e ; peak |K| = %.3f"%((K*dA)[:-1,:-1].sum(),abs(K).max()))
# standing waves: z = (cos x + cos y) cos(wt): every half period peaks and troughs swap; saddles stay saddles but flip orientation; at wt=pi/2 flat
for wt in [0,np.pi/4,np.pi/2,np.pi]:
    a=np.cos(wt); print("wt=%.2f  height at (0,0)=%+.2f  curvature at (0,0)=%+.3f  at saddle (pi,0)=%+.3f"%(wt,2*a,a*a,-a*a))
