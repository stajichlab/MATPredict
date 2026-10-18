import glob, re, sys, pickle
from Bio import SeqIO
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMMFile
K="/bigdata/stajichlab/jstajich/projects/kom"
aa=Alphabet.amino()
def fa(path):
    return list(SeqIO.parse(path,"fasta"))
# ---------- sequences
S={}  # name -> seq
def add(name,seq): S[name]=str(seq).replace("*","")
dbp={}
for f in glob.glob(K+"/db/Ascomycota/*/*/proteins.faa")+["dipo_prot.faa"]:
    tag=f.split("/")[-2] if "dipo" not in f else "dipo"
    for r in fa(f):
        m=re.search(r"name=([^|]+)",r.id)
        tg=r.id.split("|")[0] if tag=="dipo" else tag
        dbp[(tg,m.group(1) if m else r.id)]=str(r.seq)
for g in ["MATa1","MATa2","MATalpha1","MATalpha2"]:
    for (t,n),s in dbp.items():
        if t.startswith("460519_cbs") and n==g: add("Kp_"+g,s)
for (t,n),s in dbp.items():
    if t.startswith("1173061_clib918") and n=="MATA": add("Geo_MATA",s)
    if t.startswith("1173061_cbs615") and n=="MATB": add("Geo_MATB",s)
    if t.startswith("4952_w29") and n in("mata1","mata2"): add("Yl_"+n,s)
    if t.startswith("4952_clib122") and n in("matb1","matb2"): add("Yl_"+n+"_db",s)
    if t.startswith("4959_cbs767") and n.startswith("MTL"): add("Dh_"+n,s)
    if t.startswith("5476_sc5314_MTL_A") and n.startswith("MTL"): add("Ca_db_"+n,s)
    if t.startswith("5476_sc5314_MTL_alpha") and n.startswith("MTL"): add("Ca_db_"+n,s)
    if t.startswith("36911_cl143_MTL_A") and n.startswith("MTL"): add("Cl_"+n,s)
    if t.startswith("28985_nrrl-y-1140_MATsc_MATa") and n.startswith("MAT"): add("Kl_db_"+n,s)
    if t.startswith("4932_s288c_MATsc_MATalpha") and n.startswith("MAT"): add("Sc_db_"+n,s)
    if t.startswith("4932_s288c_MATsc_HMRa") and n.startswith("MAT"): add("Sc_HMR_"+n,s)
    if t.startswith("4956_da2") and n.startswith("alpha"): add("Zr_"+n,s)
    if t.startswith("381046_cbs-6340") and n.startswith("MAT"): add("Lk_"+n,s)
    if "MAT1-2-1" == n and "MAT_MAT1-2" in t and "Pezizo" not in t and "Pz" not in t:
        S.setdefault("Pez_MAT1-2-1_"+t,str(s))
pez=[k for k in S if k.startswith("Pez_")][:3]
# UniProt
up={}
for f in glob.glob("up/up_*.faa"):
    for r in fa(f):
        up[r.id.split("|")[1]]=(r.description,str(r.seq))
named={"Sc_MATALPHA1":"P0CY06","Sc_MATALPHA2":"P0CY08","Sc_HMRA1":"P0CY11","Sc_HMRA2":"P0CY13","Kl_MATA1":"Q9HG12","Kl_MATA2":"Q9HG13",
"Kl_ALPHA1":"Q08398","Kl_ALPHA2":"Q9HDS5","Ca_MTLA1":"Q9UW23","Ca_MTLA2":"Q71U11","Ca_ALPHA1":"Q9UW19","Ca_ALPHA2":"Q9UW22","Yl_MATB2_up":"Q707Y0",
"Cg_MATA1":"Q708A7","Vp_HMRa1":"A7TFG3","Pk_ALPHA1":"A0A345Z1E9","Pk_ALPHA2":"A0A345Z1F0","Op_A1_a":"W1QBK4","Op_A1_b":"W1QBL1","Op_ALPHA2":"W1QES5"}
for n,a in named.items(): add(n,up[a][1])
for a,n in [("P11632","ctl_NHP6A"),("Q05066","ctl_SRY"),("P33417","ctl_IXR1")]:
    r=fa(f"ctl/{a}.faa")[0]; add(n,r.seq)
for n in pez: pass
print("n named",len(S))
for n,s in S.items(): print("LEN",n,len(s))
pickle.dump(S,open("S.pkl","wb"))
# ---------- profiles
hmms=[]
for f in sorted(glob.glob("pf/*.hmm"))+glob.glob(K+"/db/Ascomycota/classifiers/MAT_core/*.hmm"):
    with HMMFile(f) as h:
        for x in h: hmms.append(x)
q=[TextSequence(name=n.encode(),sequence=s).digitize(aa) for n,s in S.items()]
res={}
for h in hmms:
    for hits in pyhmmer.hmmsearch(h,q,E=10):
        for hit in hits:
            d=hit.best_domain
            res[(str(hit.name),str(h.name))]="%.1f/%.0e@%d-%d"%(hit.score,hit.evalue,d.alignment.target_from,d.alignment.target_to)
pn=[str(h.name) for h in hmms]
print("PROFILES",pn)
print("PROFILE_TABLE\tprotein\t"+"\t".join(pn))
for n in S:
    print("PT\t"+n+"\t"+"\t".join(res.get((n,p),"-") for p in pn))
for h in hmms:
    print("GA",h.name,h.cutoffs.gathering, h.cutoffs.trusted)
# ---------- jackhmmer on pooled DB
pool=[]
seen=set()
for a,(d,s) in up.items():
    pool.append(TextSequence(name=a.encode(),description=d.encode(),sequence=s).digitize(aa)); seen.add(a)
for (t,n),s in dbp.items():
    pool.append(TextSequence(name=("db|"+t+"|"+n).encode(),sequence=s.replace("*","")).digitize(aa))
print("pool",len(pool))
for qn in ["Kp_MATa2","Kp_MATa1","Geo_MATA","Kl_MATA2"]:
    query=TextSequence(name=qn.encode(),sequence=S[qn]).digitize(aa)
    print("##JACK",qn)
    for it in pyhmmer.jackhmmer(query,pool,max_iterations=5,incE=1e-3,E=1e-3):
        hits=it.hits
        inc=[h for h in hits if h.included]
        print(" iter",it.iteration,"hits",len(hits),"included",len(inc),"converged",it.converged)
        last=(hits,inc)
    hits,inc=last
    for h in list(inc)[:40]:
        d=h.description.decode() if isinstance(h.description,bytes) else (h.description or "")
        print("   %s\t%.1f\t%.0e\t%s"%(str(h.name),h.score,h.evalue,d[:90]))
