"""Every-particle-to-every-particle relational data (TRPG), built 13 Nov 2026.
Values: PDG 2024/2025 review (masses, CKM magnitudes). Quark masses: MS-bar (u,d,s at 2 GeV; c,b at own mass); top = pole.
"""
import numpy as np, pandas as pd, itertools, math
P = [ # name, family, mass MeV, charge, color(1=colored), spin, weak isospin T3 (left-handed), generation
("electron","lepton",0.51099895,-1,0,0.5,-0.5,1),("muon","lepton",105.6583755,-1,0,0.5,-0.5,2),("tau","lepton",1776.93,-1,0,0.5,-0.5,3),
("nu_e","neutrino",0.0,0,0,0.5,0.5,1),("nu_mu","neutrino",0.0,0,0,0.5,0.5,2),("nu_tau","neutrino",0.0,0,0,0.5,0.5,3),
("up","quark",2.16,2/3,1,0.5,0.5,1),("down","quark",4.70,-1/3,1,0.5,-0.5,1),("strange","quark",93.5,-1/3,1,0.5,-0.5,2),
("charm","quark",1273.0,2/3,1,0.5,0.5,2),("bottom","quark",4183.0,-1/3,1,0.5,-0.5,3),("top","quark",172570.0,2/3,1,0.5,0.5,3),
("photon","gauge boson",0.0,0,0,1,0,0),("gluon","gauge boson",0.0,0,1,1,0,0),("W","gauge boson",80369.2,1,0,1,1,0),
("Z","gauge boson",91188.0,0,0,1,0,0),("Higgs","scalar",125200.0,0,0,0,-0.5,0),
("proton","baryon",938.27208816,1,0,0.5,0.5,0),("neutron","baryon",939.56542052,0,0,0.5,-0.5,0)]
cols=["name","family","mass_MeV","charge","colored","spin","T3","gen"]
df=pd.DataFrame(P,columns=cols)
CKM={("up","down"):0.97367,("up","strange"):0.22431,("up","bottom"):0.00382,("charm","down"):0.221,("charm","strange"):0.975,
     ("charm","bottom"):0.0408,("top","down"):0.0086,("top","strange"):0.0415,("top","bottom"):0.999}
lep={("electron","nu_e"),("muon","nu_mu"),("tau","nu_tau")}
v=246219.65  # Higgs vacuum value, MeV
def W_amp(a,b):
    k=(a,b) if (a,b) in CKM else (b,a)
    if k in CKM: return CKM[k]
    if (a,b) in lep or (b,a) in lep: return 1.0
    return 0.0
rows=[]
for i,j in itertools.combinations(range(len(df)),2):
    A,B=df.iloc[i],df.iloc[j]
    mi,mj=A.mass_MeV,B.mass_MeV
    r={"particle_i":A["name"],"particle_j":B["name"],"m_i":mi,"m_j":mj}
    if mi>0 and mj>0:
        r["mass_ratio"]=max(mi,mj)/min(mi,mj); r["M_ij"]=(mi+mj)/math.sqrt(mi*mj); r["Y_ij"]=math.log10(r["M_ij"]); r["S_ij"]=0.5*math.log(mi/mj)
    qq=A.charge*B.charge
    r["charge_product"]=round(qq,4)
    r["EM_push_pull"]="push (repel)" if qq>0 else ("pull (attract)" if qq<0 else "none")
    r["strong_link"]=int(A.colored and B.colored)
    r["W_transition_amplitude"]=W_amp(A["name"],B["name"])
    r["T3_product"]=A.T3*B.T3
    yi=math.sqrt(2)*mi/v if A.family not in("gauge boson",) else 0; yj=math.sqrt(2)*mj/v if B.family not in("gauge boson",) else 0
    r["Higgs_coupling_product"]=yi*yj
    r["same_generation"]=int(A.gen==B.gen and A.gen>0)
    rows.append(r)
pairs=pd.DataFrame(rows)
# Test: Cabibbo-type relations |V| vs sqrt(mass ratio) = exp(S)
tests=[]
for (u,d),V in CKM.items():
    if u=="up" and d=="down" or u=="charm" and d=="strange" or u=="top" and d=="bottom": continue
    m=dict(zip(df.name,df.mass_MeV))
    cands={"sqrt(m_lighter_down / m_heavier_down)":None}
    tests.append((u,d,V))
gst=[("V_us","sqrt(m_d/m_s)",CKM[("up","strange")],math.sqrt(4.70/93.5)),
     ("V_ub","sqrt(m_u/m_t)",CKM[("up","bottom")],math.sqrt(2.16/172570)),
     ("V_cb","sqrt(m_s/m_b)",CKM[("charm","bottom")],math.sqrt(93.5/4183)),
     ("V_cb","sqrt(m_c/m_t)",CKM[("charm","bottom")],math.sqrt(1273/172570)),
     ("V_td","sqrt(m_d/m_b)",CKM[("top","down")],math.sqrt(4.70/4183)),
     ("V_ts","sqrt(m_s/m_b)",CKM[("top","strange")],math.sqrt(93.5/4183))]
gt=pd.DataFrame(gst,columns=["CKM element","relational guess = exp(S_ij)","measured","predicted"])
gt["percent_off"]=(gt.predicted/gt.measured-1)*100
print(gt.round(4).to_string(index=False))
with pd.ExcelWriter("TRPG_Every_Particle_to_Every_Particle.xlsx") as w:
    df.to_excel(w,sheet_name="Particles",index=False); pairs.to_excel(w,sheet_name="All_Pairs",index=False); gt.to_excel(w,sheet_name="Mixing_vs_Mass_Test",index=False)
pairs.to_csv("TRPG_All_Pairs.csv",index=False)
print(len(pairs),"pairs")
# summary: how many pairs linked by each force
print("EM-linked:",(pairs.charge_product!=0).sum()," strong-linked:",pairs.strong_link.sum()," W-linked:",(pairs.W_transition_amplitude>0).sum())

# ---- Scale_Views (added 13 Nov 2026, the idea: each particle sets its own scale view of the one universe) ----
# Constants: CODATA 2018; H0 = 67.4 km/s/Mpc (Planck 2018).
hbar=1.054571817e-34;c=2.99792458e8;G=6.67430e-11;H0=67.4e3/3.0857e22;MeVJ=1.602176634e-13
mPl=math.sqrt(hbar*c/G); RH=c/H0; Muni=c**3/(2*G*H0)
massive=df[df.mass_MeV>0].reset_index(drop=True)
view=pd.DataFrame({"reference_particle":massive.name})
for _,B in massive.iterrows(): view["%s_in_ref_units"%B["name"]]=massive.mass_MeV.map(lambda mr:B.mass_MeV/mr)
kg=massive.mass_MeV*MeVJ/c**2
view.insert(1,"Compton_length_fm",hbar/(kg*c)*1e15)
view.insert(2,"gravity_vs_quantum_1_over_alphaG",(mPl/kg)**2)
view.insert(3,"Hubble_radius_in_ref_lengths",RH/(hbar/(kg*c)))
view.insert(4,"universe_mass_in_ref_particles",Muni/kg)
with pd.ExcelWriter("TRPG_Every_Particle_to_Every_Particle.xlsx",mode="a",engine="openpyxl",if_sheet_exists="replace") as w:
    view.to_excel(w,sheet_name="Scale_Views",index=False)
print(view.iloc[:,:5].to_string())
