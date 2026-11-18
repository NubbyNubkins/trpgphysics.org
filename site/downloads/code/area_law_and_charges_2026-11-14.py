"""Finishing calculations (14 Nov 2026).
A) Area law: entanglement of a region of a wave medium (2D lattice of coupled oscillators, ground state) grows with its BOUNDARY, not its volume
   (Bombelli, Koul, Lee & Sorkin 1986; Srednicki 1993). This is the input Jacobson (1995) needs to derive Einstein's equations.
B) Charges from consistency: anomaly cancellation fixes the Standard Model charges and the number of colours."""
import numpy as np, math
N=36; m=0.01
idx=lambda i,j:i*N+j
K=np.zeros((N*N,N*N))
for i in range(N):
    for j in range(N):
        a=idx(i,j); K[a,a]=4+m*m
        for di,dj in((1,0),(0,1)):
            ii,jj=i+di,j+dj
            if ii<N and jj<N: b=idx(ii,jj); K[a,b]=K[b,a]=-1
w,U=np.linalg.eigh(K)
X=0.5*(U/np.sqrt(w))@U.T; P=0.5*(U*np.sqrt(w))@U.T
res=[]
for L in range(2,15,2):
    o=(N-L)//2; A=[idx(i,j) for i in range(o,o+L) for j in range(o,o+L)]
    ev=np.linalg.eigvals(X[np.ix_(A,A)]@P[np.ix_(A,A)]).real; nu=np.sqrt(np.clip(ev,0.25+1e-12,None))
    S=sum((n+.5)*math.log(n+.5)-(n-.5)*math.log(n-.5) for n in nu if n>0.5+1e-9)
    res.append((L,4*L,L*L,S)); print("block %2dx%-2d boundary %3d volume %4d  entanglement %.3f  per boundary %.4f"%(L,L,4*L,L*L,S,S/(4*L)))
L=np.array([r[0] for r in res]); S=np.array([r[3] for r in res])
cb=np.polyfit(4*L,S,1); cv=np.polyfit(L*L,S,1)
rb=np.sqrt(np.mean((np.polyval(cb,4*L)-S)**2)); rv=np.sqrt(np.mean((np.polyval(cv,L*L)-S)**2))
print("fit to boundary: error %.3f ; fit to volume: error %.3f -> %s law"%(rb,rv,"AREA (boundary)" if rb<rv else "volume"))
print("entanglement per unit boundary ~ %.3f -> Jacobson/Bekenstein: G = a^2/(4 x that) in medium units: G ~ %.2f a^2 (gravity strength set by the medium's spacing a)"%(cb[0],1/(4*cb[0])))
# B) charges
print("\nB) Charges from consistency (anomaly cancellation; Geng & Marshak 1989; Minahan, Ramond & Warner 1990; Bar & Wiese 2001).")
print("   One generation with N colours: quark doublet Q(N,2), u(N), d(N), lepton doublet L(2), e(1); fix the scale with Y_L = -1/2.")
import sympy as sp
YQ, Yu, Yd, Ye, N = sp.symbols("Y_Q Y_u Y_d Y_e N")
YL = sp.Rational(-1, 2)
for Nc in (1, 2, 3, 4, 5):
    eqs = [2*YQ - Yu - Yd,                                   # SU(N)^2 U(1)
           Nc*YQ + YL,                                       # SU(2)^2 U(1)
           Nc*(2*YQ - Yu - Yd) + 2*YL - Ye,                  # gravity^2 U(1)
           Nc*(2*YQ**3 - Yu**3 - Yd**3) + 2*YL**3 - Ye**3]   # U(1)^3
    sols = sp.solve(eqs, [YQ, Yu, Yd, Ye], dict=True)
    for so in sols:
        qu = sp.Rational(1, 2) + so[YQ]; qd = -sp.Rational(1, 2) + so[YQ]
        if so[Yu] > so[Yd]:
            print(f"   N = {Nc}: up-quark charge {qu}, down-quark charge {qd}, electron {so[Ye] + 0}, neutrino 0")
print("   Result: for ANY number of colours N the anomalies cancel, with up = 1/2 + 1/(2N) and down = -1/2 + 1/(2N).")
print("   So the charges are forced ONCE N is known, and N = 3 gives exactly +2/3 and -1/3. N itself is not forced by")
print("   anomalies; it is measured (e.g. the neutral pion decay rate and the e+e- -> hadrons ratio both need N = 3).")
