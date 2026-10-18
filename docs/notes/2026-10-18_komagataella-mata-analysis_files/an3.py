import glob, pickle, subprocess, collections
from Bio import SeqIO
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, TextMSA
from pyhmmer.plan7 import HMMFile, Builder, Background
K="/bigdata/stajichlab/jstajich/projects/kom"
aa=Alphabet.amino()
S=pickle.load(open("S.pkl","rb"))
up={}
for f in glob.glob("up/up_*.faa"):
    for r in SeqIO.parse(f,"fasta"): up[r.id.split("|")[1]]=(r.description,str(r.seq))
desc={a:d.split(" OS=")[0][:70] for a,(d,s) in up.items()}
ext={"ext|Dh_MTLA2":S["Dh_MTLA2"],"ext|Cl_MTLA2":S["Cl_MTLA2"],"ext|Geo_MATA":S["Geo_MATA"],"ext|Yl_mata2":S["Yl_mata2"],
     "ext|Cl_MTLA1":S["Cl_MTLA1"],"ext|Dh_MTLA1":S["Dh_MTLA1"],"ext|Yl_mata1":S["Yl_mata1"],"ext|Geo_MATB":S["Geo_MATB"]}
pool=[TextSequence(name=a.encode(),sequence=s).digitize(aa) for a,(d,s) in up.items()]+[TextSequence(name=n.encode(),sequence=s).digitize(aa) for n,s in ext.items()]
# which UniProt accession = which named protein
acc={"Kp_MATa2":"A0A1G4KQN4","Kl_MATA2":"Q9HG13","Ca_MTLA2":"Q71U11","Kp_MATa1":"F2QXQ4","Kl_MATA1":"Q9HG12","Ca_MTLA1":"Q9UW23","Sc_HMRA1":"P0CY11","Cg_MATA1":"Q708A7","Vp_HMRa1":"A7TFG3","Op_A1_a":"W1QBK4"}
extn={"Dh_MTLA2":"ext|Dh_MTLA2","Cl_MTLA2":"ext|Cl_MTLA2","Geo_MATA":"ext|Geo_MATA","Yl_mata2":"ext|Yl_mata2","Cl_MTLA1":"ext|Cl_MTLA1","Dh_MTLA1":"ext|Dh_MTLA1","Yl_mata1":"ext|Yl_mata1"}
poolname=lambda n: acc.get(n) or extn[n]
def build(names,tag):
    with open("t.faa","w") as f:
        for n in names: f.write(f">{n}\n{S[n]}\n")
    r=subprocess.run(["mafft","--localpair","--maxiterate","1000","--quiet","t.faa"],capture_output=True,text=True)
    aln=collections.OrderedDict();cur=None
    for l in r.stdout.splitlines():
        if l.startswith(">"): cur=l[1:].strip(); aln[cur]=""
        else: aln[cur]+=l.strip()
    open(f"{tag}.aln.faa","w").write("".join(f">{n}\n{s}\n" for n,s in aln.items()))
    msa=TextMSA(name=tag.encode(),sequences=[TextSequence(name=n.encode(),sequence=s) for n,s in aln.items()]).digitize(aa)
    b=Builder(aa); bg=Background(aa)
    hmm,_,_=b.build_msa(msa,bg)
    return hmm
def search(hmm,posnames):
    hits=list(pyhmmer.hmmsearch(hmm,pool,E=1000))[0]
    ranked=[(str(h.name),h.score,h.evalue) for h in hits]
    return ranked
def report(title,groups,negatives_note):
    print("=====",title)
    allnames=[n for g in groups.values() for n in g]
    full=build(allnames,title.split()[0]+"_all")
    pos={poolname(n) for n in allnames}
    r=search(full,pos)
    print(" ALL-lineage profile (%d seqs): top hits"%len(allnames))
    for n,s,e in r[:14]: print("    %-14s %6.1f %.0e %s %s"%(n,s,e,"POS" if n in pos else "   ",desc.get(n,"")))
    nonpos=[x for x in r if x[0] not in pos]
    print("    best non-positive: %.1f (%s)"%(nonpos[0][1],desc.get(nonpos[0][0],nonpos[0][0])))
    for gname,gm in groups.items():
        train=[n for n in allnames if n not in gm]
        h=build(train,title.split()[0]+"_LO_"+gname)
        r=search(h,pos)
        rank={n:i for i,(n,s,e) in enumerate(r)}
        pn={poolname(m) for m in gm}
        npos=[(n,s,e) for n,s,e in r if n not in {poolname(x) for x in allnames}]
        print(" LEAVE OUT %-12s trained on %d seqs"%(gname,len(train)))
        for m in gm:
            pm=poolname(m)
            if pm in rank:
                n,s,e=r[rank[pm]]
                nneg=sum(1 for x in r[:rank[pm]] if x[0] not in {poolname(y) for y in allnames})
                print("     held-out %-10s score %6.1f E %.0e  rank %d (non-positives above: %d)"%(m,s,e,rank[pm]+1,nneg))
            else: print("     held-out %-10s no hit at E<=1000"%m)
        print("       best non-positive above/at E<=1000: %.1f (%s); 5th: %.1f"%(npos[0][1],desc.get(npos[0][0],npos[0][0]),npos[4][1]))
a2={"Komagataella":["Kp_MATa2"],"Saccharomycetaceae":["Kl_MATA2"],"Serinales":["Ca_MTLA2","Dh_MTLA2","Cl_MTLA2"],"Dipodascales":["Geo_MATA","Yl_mata2"]}
report("MATa2 HMG-like set",a2,"")
a1={"Komagataella":["Kp_MATa1"],"Saccharomycetaceae":["Kl_MATA1","Sc_HMRA1","Cg_MATA1","Vp_HMRa1"],"Serinales":["Ca_MTLA1","Dh_MTLA1","Cl_MTLA1"],"Ogataea":["Op_A1_a"],"Dipodascales":["Yl_mata1"]}
report("MATa1 homeodomain set",a1,"")
# reference profiles over held-out/all positives for comparison
print("===== comparison profiles over the same proteins (rank in pool of 72k)")
for f in sorted(glob.glob("pf/PF00505.hmm"))+sorted(glob.glob("pf/PF09011.hmm"))+sorted(glob.glob("pf/PF00046.hmm"))+glob.glob(K+"/db/Ascomycota/classifiers/MAT_core/MAT1-2-1.hmm"):
    with HMMFile(f) as hf: hmm=next(iter(hf))
    r=list(pyhmmer.hmmsearch(hmm,pool,E=1000))[0]
    rk={str(h.name):(i+1,h.score,h.evalue) for i,h in enumerate(r)}
    print(" profile",hmm.name.decode() if isinstance(hmm.name,bytes) else hmm.name,"n hits E<=1000:",len(rk))
    for n in ["Kp_MATa2","Kl_MATA2","Ca_MTLA2","Dh_MTLA2","Cl_MTLA2","Geo_MATA","Yl_mata2","Kp_MATa1","Kl_MATA1","Yl_mata1"]:
        pm=poolname(n)
        print("     %-10s %s"%(n,("rank %d score %.1f E %.0e"%rk[pm]) if pm in rk else "no hit E<=1000"))
