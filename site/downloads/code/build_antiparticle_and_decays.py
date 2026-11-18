"""Antiparticle matrices and decay data (14 Nov 2026), added to TRPG_Every_Particle_to_Every_Particle.xlsx.
Sheets: Antiparticles, AntiAnti_Pairs, Particle_vs_Antiparticle (19x19 = 361), Decays (PDG 2024), Decay_Relations (Fermi-theory check).
Annihilation rule: a pair can annihilate into neutral force carriers if every conserved number sums to zero
(charge, baryon number, each lepton-family number, quark flavour; colour must be able to cancel, so a coloured and a
colourless particle cannot annihilate). Proton-antineutron pairs annihilate into charged pions. It can fuse into a W if the pair
has charge +/-1, zero baryon number and is a weak-doublet partner pair (CKM or lepton-family); otherwise only scattering."""
import pandas as pd, itertools, math
F="TRPG_Every_Particle_to_Every_Particle.xlsx"
P=pd.read_excel(F,sheet_name="Particles")
B={"up":1/3,"down":1/3,"strange":1/3,"charm":1/3,"bottom":1/3,"top":1/3,"proton":1,"neutron":1}
Lf={"electron":"e","nu_e":"e","muon":"mu","nu_mu":"mu","tau":"tau","nu_tau":"tau"}
self_conj={"photon","Z","Higgs","gluon"}   # gluon octet as a set is self-conjugate
CKM={("up","down"):0.97367,("up","strange"):0.22431,("up","bottom"):0.00382,("charm","down"):0.221,("charm","strange"):0.975,
     ("charm","bottom"):0.0408,("top","down"):0.0086,("top","strange"):0.0415,("top","bottom"):0.999}
def qn(name,anti):
    s=-1 if anti and name not in self_conj else 1
    r=P.set_index("name").loc[name]
    L={"e":0,"mu":0,"tau":0}
    if name in Lf: L[Lf[name]]=1*s
    return dict(Q=r.charge*s,B=B.get(name,0)*s,**{"L_"+k:v for k,v in L.items()})
def label(n,anti): return ("anti-"+n if anti and n not in self_conj else n)
anti=P.copy(); anti["name"]=[label(n,True) for n in P.name]; anti["charge"]=[-c if n not in self_conj else c for n,c in zip(P.name,P.charge)]
anti["T3"]=[-t for t in P.T3]; anti["baryon_number"]=[-B.get(n,0) for n in P.name]
def pair_row(a,aa,b,ba):
    A,Bq=qn(a,aa),qn(b,ba); tot={k:A[k]+Bq[k] for k in A}
    ma=P.set_index("name").loc[a].mass_MeV; mb=P.set_index("name").loc[b].mass_MeV
    direct = (a==b and aa!=ba) or (a==b and a in self_conj)
    quarks={"up","down","strange","charm","bottom","top"}; nucleons={"proton","neutron"}
    col=lambda n: bool(P.set_index("name").loc[n].colored)
    if col(a)!=col(b):
        fate="no annihilation; scatter only (colour cannot cancel)"
    elif all(abs(v)<1e-9 for v in tot.values()) and a in quarks and b in quarks and a!=b:
        fate="flavour-changing: no direct annihilation (only rare loop processes)"
    elif all(abs(v)<1e-9 for v in tot.values()):
        fate="can annihilate -> neutral carriers (photons, gluons, Z, Higgs)"
    elif a in nucleons and b in nucleons and abs(tot["B"])<1e-9:
        fate="can annihilate -> pions (strong force; hadronic)"
    elif abs(abs(tot["Q"])-1)<1e-9 and abs(tot["B"])<1e-9 and (((a,b) in CKM or (b,a) in CKM) or (Lf.get(a)==Lf.get(b) and a in Lf)) and all(abs(tot[k])<1e-9 for k in ("L_e","L_mu","L_tau")):
        fate="can fuse -> W boson (weak charged current)"
    else: fate="no annihilation; scatter only"
    qq=A["Q"]*Bq["Q"]
    return {"i":label(a,aa),"j":label(b,ba),"direct_opposites":int(direct),"net_charge":round(tot["Q"],4),"net_baryon":round(tot["B"],4),
            "net_lepton_e":tot["L_e"],"net_lepton_mu":tot["L_mu"],"net_lepton_tau":tot["L_tau"],
            "EM_push_pull":"push (repel)" if qq>0 else ("pull (attract)" if qq<0 else "none"),
            "min_energy_to_make_pair_MeV":ma+mb,"fate":fate}
names=list(P.name)
aa=pd.DataFrame([pair_row(a,True,b,True) for a,b in itertools.combinations(names,2)])
pa=pd.DataFrame([pair_row(a,False,b,True) for a in names for b in names])
# Decays (PDG 2024 values)
dec=[("muon","e nu_mu anti-nu_e",1.0,"lifetime 2.1969811e-6 s"),("tau","mu nu nu",0.1739,"lifetime 290.3e-15 s"),("tau","e nu nu",0.1782,""),("tau","hadrons nu",0.648,""),
     ("W","e nu",0.1071,"width 2.085 GeV"),("W","mu nu",0.1063,""),("W","tau nu",0.1138,""),("W","hadrons (quark pairs)",0.6741,""),
     ("Z","e e",0.03363,"width 2.4955 GeV"),("Z","mu mu",0.03366,""),("Z","tau tau",0.03370,""),("Z","invisible (neutrino pairs)",0.2000,""),("Z","hadrons",0.6991,""),
     ("Higgs","b anti-b",0.582,"SM width ~4.1 MeV"),("Higgs","W W*",0.214,""),("Higgs","gluon gluon",0.082,""),("Higgs","tau tau",0.0627,""),("Higgs","c anti-c",0.0289,""),("Higgs","Z Z*",0.0262,""),("Higgs","photon photon",0.00227,""),
     ("top","W b",0.998,"width ~1.42 GeV"),("neutron","p e anti-nu_e",1.0,"lifetime 878.4 s"),("proton","(stable; > 1e34 yr)",0.0,"")]
decays=pd.DataFrame(dec,columns=["parent","products","branching_fraction","total_lifetime_or_width"])
# Fermi theory check: muon lifetime and tau lifetime from the muon via mass ratio^5 = e^(10 S)
GF=1.1663788e-5*1e-6; hbar=6.582119569e-22  # GF in MeV^-2, hbar in MeV s
mmu,mtau=105.6583755,1776.93
tau_mu=hbar/(GF**2*mmu**5/(192*math.pi**3)); S=0.5*math.log(mtau/mmu)
tau_tau=tau_mu*0.1782*math.exp(-10*S)
rel=pd.DataFrame([{"quantity":"muon lifetime from Fermi theory","value":tau_mu,"measured":2.1969811e-6},
                  {"quantity":"tau lifetime = muon lifetime x BR(tau->e nu nu) x e^(-10 S_tau,mu)","value":tau_tau,"measured":290.3e-15}])
with pd.ExcelWriter(F,mode="a",engine="openpyxl",if_sheet_exists="replace") as w:
    anti.to_excel(w,sheet_name="Antiparticles",index=False); aa.to_excel(w,sheet_name="AntiAnti_Pairs",index=False)
    pa.to_excel(w,sheet_name="Particle_vs_Antiparticle",index=False); decays.to_excel(w,sheet_name="Decays",index=False); rel.to_excel(w,sheet_name="Decay_Relations",index=False)
print(pa.fate.value_counts().to_string()); print("direct opposites:",pa.direct_opposites.sum())
print(rel.to_string(index=False))
print(pa[pa.fate.str.contains("W")][["i","j"]].to_string(index=False))
