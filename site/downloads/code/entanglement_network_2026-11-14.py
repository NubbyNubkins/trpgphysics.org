"""M15: 'degrees of separation' for particles (14 Nov 2026).
1) Small-world estimate L ~ ln N / ln k (Watts & Strogatz 1998) for people and for the universe's particles.
2) Degrees of separation in the every-particle table (K5): links = any direct interaction."""
import math, pandas as pd, itertools
for name,N,k in [("people (8e9, ~150 close contacts; Dunbar)",8e9,150),("particles in visible universe (~1e80), 10 links each",1e80,10),("same, 1000 links each",1e80,1000)]:
    print("%-58s ~%.1f steps"%(name,math.log(N)/math.log(k)))
df=pd.read_excel("TRPG_Every_Particle_to_Every_Particle.xlsx",sheet_name="All_Pairs")
def graph(use_higgs):
    E={}
    for _,r in df.iterrows():
        linked = r.charge_product!=0 or r.strong_link==1 or r.W_transition_amplitude>0 or (use_higgs and r.Higgs_coupling_product>0)
        if linked: E.setdefault(r.particle_i,set()).add(r.particle_j); E.setdefault(r.particle_j,set()).add(r.particle_i)
    return E
names=sorted(set(df.particle_i)|set(df.particle_j))
for use_higgs in (False,True):
    E=graph(use_higgs); dist=[]; unreached=0
    for a in names:
        seen={a:0}; fr=[a]
        while fr:
            nx=[]
            for u in fr:
                for v in E.get(u,()):
                    if v not in seen: seen[v]=seen[u]+1; nx.append(v)
            fr=nx
        for b in names:
            if b>a: (dist.append(seen[b]) if b in seen else None)
            if b>a and b not in seen: unreached+=1
    print("K5 table, Higgs links %s: max steps %d, average %.2f, pairs never connected %d"%("on" if use_higgs else "off",max(dist),sum(dist)/len(dist),unreached))
# Add force-carrier couplings (photon<->charged, gluon<->quarks, W/Z<->all weak-charged, Higgs<->massive)
P=pd.read_excel("TRPG_Every_Particle_to_Every_Particle.xlsx",sheet_name="Particles").set_index("name")
E=graph(True)
def link(a,b): E.setdefault(a,set()).add(b); E.setdefault(b,set()).add(a)
for n,r in P.iterrows():
    if n in("photon","gluon","W","Z","Higgs"): continue
    if r.charge!=0: link("photon",n)
    if r.colored: link("gluon",n)
    if r.T3!=0 or r.charge!=0: link("Z",n); link("W",n)
    if r.mass_MeV>0: link("Higgs",n)
link("W","photon"); link("W","Z"); link("Higgs","W"); link("Higgs","Z")
dist=[]; un=0
for a in names:
    seen={a:0}; fr=[a]
    while fr:
        nx=[]
        for u in fr:
            for v in E.get(u,()):
                if v not in seen: seen[v]=seen[u]+1; nx.append(v)
        fr=nx
    for b in names:
        if b>a:
            if b in seen: dist.append(seen[b])
            else: un+=1
print("with force carriers: max steps %d, average %.2f, never connected %d (of 171)"%(max(dist),sum(dist)/len(dist),un))
