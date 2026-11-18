"""Quark masses from linked Koide triples ('waterfall', Rivero 2011; cf. Rodejohann & Zhang 2011). 14 Nov 2026.
Each overlapping triple (with one negative sqrt where needed) is set to Q=2/3; start from top and bottom only.
Caveat: quark masses depend on scheme/scale (MS-bar), so % comparisons are approximate."""
import math
from scipy.optimize import brentq
m={"u":2.16,"d":4.70,"s":93.5,"c":1273,"b":4183,"t":172570}
def Q(vals,signs): r=[s*math.sqrt(v) for v,s in zip(vals,signs)]; return sum(vals)/sum(r)**2
print("Q(c,b,t)=%.4f  Q(-s,c,b)=%.4f  Q(d,s,-c)? see chain"%(Q([m['c'],m['b'],m['t']],[1,1,1]),Q([m['s'],m['c'],m['b']],[-1,1,1])))
def solve_light(heavy2,heavy1,sign_light):
    f=lambda x: Q([x,heavy2,heavy1],[sign_light,1,1])-2/3
    xs=[10**(k/200) for k in range(-600,1200)]; out=[]
    for a,b in zip(xs,xs[1:]):
        if f(a)*f(b)<0: out.append(brentq(f,a,b))
    return out
c=solve_light(m['b'],m['t'],1); print("charm from (b,t):",[round(v,1) for v in c],"measured",m['c'])
cc=min(c,key=lambda v:abs(v-m['c']))
s=solve_light(cc,m['b'],-1); print("strange from (c,b), negative root:",[round(v,1) for v in s],"measured",m['s'])
