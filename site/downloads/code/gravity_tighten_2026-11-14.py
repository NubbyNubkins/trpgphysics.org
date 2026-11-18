"""Tightening the gravity chain (14 Nov 2026).
A) Area law in the crossing-wave medium: couplings vary as 1 + 0.5 cos(x) cos(y) (resistance from crossing waves).
B) Equal action on all layers: the area-law coefficient for fields of very different mass (different 'layers')."""
import numpy as np, math
def ent(N,m,crossing,L):
    idx=lambda i,j:i*N+j; K=np.zeros((N*N,N*N))
    for i in range(N):
        for j in range(N):
            a=idx(i,j)
            for di,dj in((1,0),(0,1)):
                ii,jj=i+di,j+dj
                if ii<N and jj<N:
                    xm,ym=(i+ii)/2,(j+jj)/2
                    k=1+0.5*math.cos(2*math.pi*xm/6)*math.cos(2*math.pi*ym/6) if crossing else 1.0
                    b=idx(ii,jj); K[a,b]=K[b,a]=-k; K[a,a]+=k; K[b,b]+=k
            K[a,a]+=m*m
    w,U=np.linalg.eigh(K); X=0.5*(U/np.sqrt(w))@U.T; P=0.5*(U*np.sqrt(w))@U.T
    out=[]
    for l in L:
        o=(N-l)//2; A=[idx(i,j) for i in range(o,o+l) for j in range(o,o+l)]
        nu=np.sqrt(np.clip(np.linalg.eigvals(X[np.ix_(A,A)]@P[np.ix_(A,A)]).real,0.25+1e-12,None))
        out.append(sum((n+.5)*math.log(n+.5)-(n-.5)*math.log(n-.5) for n in nu if n>0.5+1e-9))
    return np.array(out)
N=36; L=np.array([4,6,8,10,12,14])
S=ent(N,0.01,True,L)
eb=np.sqrt(np.mean((np.polyval(np.polyfit(4*L,S,1),4*L)-S)**2)); ev=np.sqrt(np.mean((np.polyval(np.polyfit(L*L,S,1),L*L)-S)**2))
print("A) crossing-wave medium: S =",np.round(S,3)," boundary-fit error %.3f vs volume-fit error %.3f (%.0fx better)"%(eb,ev,ev/eb))
print("B) area coefficient (entanglement per boundary unit) for different layer masses:")
for m in (0.01,0.1,0.3,1.0):
    S=ent(N,m,False,L); c=np.polyfit(4*L,S,1)[0]; print("   mass %.2f: coefficient %.4f"%(m,c))
