"""The Higgs mass from the grain of space: the asymptotic-safety rule tested inside the framework.
Rule from Shaposhnikov & Wetterich (2010, PLB 683, 196): gravity drives the Higgs self-hold lambda to zero at the grain.
TEST FIXED BEFORE RUNNING:
  find the Higgs mass m_H for which lambda(grain) = 0, two-loop running (Machacek & Vaughn 1983-85; Luo, Wang & Xiao 2003),
  starting values at mu = m_t from the NNLO fit of Buttazzo, Degrassi, Giardino, Giudice, Sala, Salvio & Strumia (2013).
  main grain = Planck/sqrt2 = 8.63e18 GeV (the framework's grain, gives alpha_s 0.1178); also shown: Planck 1.22e19.
  runs: (a) Standard Model only; (b) with the framework's far side: 3 far-side weak pairs (spin 1/2, colourless, both-handed, Y=1/2) at 2.4e8 GeV
        and 3 D partners (3,-1/3) at the meeting F = 4.6e13 GeV; no link between the far-side pairs and the Higgs (the framework has none).
  uncertainty: top pole mass 172.57 +- 0.29 (PDG 2024) plus 0.5 GeV for the pole/Monte-Carlo mass ambiguity, alpha_s 0.1180 +- 0.0009,
  plus 0.3 GeV for the missing three-loop terms (Buttazzo et al.). Combined in quadrature.
  PASS line: within 2 sigma of the measured 125.20 +- 0.11 GeV. FAIL: beyond 3 sigma. Between: PARTIAL.
"""
import math, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
K=16*math.pi**2
MPL=1.220890e19; GRAIN=MPL/math.sqrt(2); MMG=2.4e8; F=4.6e13
MH_MEAS=125.20; dMH_MEAS=0.11
def start(mh,mt=172.57,a3=0.1180,mw=80.3692):
    lam=0.12604+0.00206*(mh-125.15)-0.00004*(mt-173.34)
    yt=0.93690+0.00556*(mt-173.34)-0.00042*(a3-0.1184)/0.0007
    g3=1.1666+0.00314*(a3-0.1184)/0.0007-0.00046*(mt-173.34)
    g2=0.64779+0.00004*(mt-173.34)+0.00011*(mw-80.384)/0.014
    g1=0.35830+0.00011*(mt-173.34)-0.00020*(mw-80.384)/0.014     # g' (SM normalisation)
    return np.array([g1,g2,g3,yt,lam])
# far-side pieces, SM normalisation of g': (db1, db2, db3) and two-loop matrices, per Dirac multiplet
dbL=np.array([2/3,2/3,0]); dBL=np.array([[25/9*9/50,5/3*9/10,0],[5/3*3/10,49/6,0],[0,0,0]])
dbD=np.array([4/9,0,2/3]); dBD=np.array([[25/9*4/75,0,5/3*16/15],[0,0,0],[5/3*2/15,0,38/3]])
b0=np.array([41/6,-19/6,-7.0]); B0=np.array([[199/18,9/2,44/3],[3/2,35/6,12],[11/6,9/2,-26]]); cT=np.array([17/6,3/2,2])
def beta(t,y,far):
    g1,g2,g3,yt,l=y; mu=172.57*math.exp(t)
    b=b0.copy(); B=B0.copy()
    if far:
        if mu>MMG: b=b+3*dbL; B=B+3*dBL
        if mu>F:   b=b+3*dbD; B=B+3*dBD
    g=np.array([g1,g2,g3]); G=g*g
    dg=g**3*b/K+g**3*(B@G-cT*yt*yt)/K**2
    a1,a2,a3,T=G[0],G[1],G[2],yt*yt
    dyt=yt*(4.5*T-8*a3-2.25*a2-17/12*a1)/K + yt*(-12*T*T+T*(36*a3+225/16*a2+131/16*a1-12*l)+6*l*l
        -108*a3*a3+9*a2*a3+19/9*a1*a3-23/4*a2*a2-3/4*a1*a2+1187/216*a1*a1)/K**2
    dl=(24*l*l-6*T*T+3/8*(2*a2*a2+(a2+a1)**2)+l*(12*T-9*a2-3*a1))/K
    dl+=(-312*l**3+36*l*l*(3*a2+a1)-144*l*l*T-3*l*T*T+80*l*a3*T+45/2*l*a2*T+85/6*l*a1*T
        -73/8*l*a2*a2+39/4*l*a1*a2+629/24*l*a1*a1+30*T**3-32*a3*T*T-8/3*a1*T*T-9/4*a2*a2*T
        +21/2*a1*a2*T-19/4*a1*a1*T+305/16*a2**3-289/48*a1*a2*a2-559/48*a1*a1*a2-379/48*a1**3)/K**2
    return [dg[0],dg[1],dg[2],dyt,dl]
def lam_at(mu,mh,far=False,**kw):
    tend=math.log(mu/172.57)
    edges=[0.0]+([math.log(MMG/172.57),math.log(F/172.57)] if far else [])+[tend]
    edges=sorted(e for e in edges if e<=tend)
    y=start(mh,**kw)
    for a,b in zip(edges[:-1],edges[1:]):
        y=solve_ivp(beta,[a,b],y,args=(far,),rtol=1e-10,atol=1e-13).y[:,-1]
    return y
def mh_needed(grain=GRAIN,far=False,**kw):
    return brentq(lambda m: lam_at(grain,m,far,**kw)[4],100,160,xtol=1e-5)
print("check against Buttazzo et al. 2013 (m_t 173.34, m_H 125.15, alpha_s 0.1184): lambda(M_Planck) =",
      round(lam_at(MPL,125.15,mt=173.34,a3=0.1184,mw=80.384)[4],4)," (their value -0.0143 at 3 loops)")
print("  and m_H for lambda(M_Planck)=0 at m_t 173.34:",round(mh_needed(MPL,mt=173.34,a3=0.1184,mw=80.384),2),
      "(their stability bound 129.6 is a nearby but not identical condition)")
print()
print("today's inputs: m_t 172.57, alpha_s 0.1180, m_W 80.369; lambda today at the grain with the measured m_H:")
for far in (False,True):
    y=lam_at(GRAIN,MH_MEAS,far); print(f"  {'with the far side     ' if far else 'Standard Model only':22s} lambda(grain) = {y[4]:+.4f}")
print()
def sigma(grain,far):
    c=mh_needed(grain,far)
    dt=(mh_needed(grain,far,mt=172.57+0.29)-mh_needed(grain,far,mt=172.57-0.29))/2
    da=(mh_needed(grain,far,a3=0.1189)-mh_needed(grain,far,a3=0.1171))/2
    dpole=abs(dt)/0.29*0.5
    s=math.sqrt(dt**2+da**2+dpole**2+0.3**2+dMH_MEAS**2)
    return c,s,dt,da,dpole
for gl,g in (("grain Planck/sqrt2 (main)",GRAIN),("grain Planck",MPL)):
    for far in (False,True):
        c,s,dt,da,dp=sigma(g,far); z=(c-MH_MEAS)/s
        verdict="PASS" if abs(z)<2 else ("FAIL" if abs(z)>3 else "PARTIAL")
        print(f"{gl:26s} {'with far side' if far else 'SM only':14s}: m_H = {c:7.2f} +- {s:.2f} GeV  (top {dt:.2f}, alpha_s {da:.2f}, pole {dp:.2f}, loops 0.30)"
              f"  vs 125.20: {z:+.1f} sigma -> {verdict}")
mt_fit=brentq(lambda mt: mh_needed(GRAIN,True,mt=mt)-MH_MEAS,165,180,xtol=1e-4)
mt_fit0=brentq(lambda mt: mh_needed(GRAIN,False,mt=mt)-MH_MEAS,165,180,xtol=1e-4)
print(f"\ntop pole mass the rule would need for m_H = 125.20: SM only {mt_fit0:.2f} GeV, with far side {mt_fit:.2f} GeV (measured 172.57 +- 0.29)")
