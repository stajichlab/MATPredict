import glob, re, sys, pickle, subprocess, os, gzip, tempfile, collections
from Bio import SeqIO
from Bio.Seq import Seq
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, TextMSA
from pyhmmer.plan7 import HMMFile, Builder, Background
K="/bigdata/stajichlab/jstajich/projects/kom"
aa=Alphabet.amino()
S=pickle.load(open("S.pkl","rb"))
up={}
for f in glob.glob("up/up_*.faa"):
    for r in SeqIO.parse(f,"fasta"):
        up[r.id.split("|")[1]]=(r.description,str(r.seq))
pool=[TextSequence(name=a.encode(),description=d.encode(),sequence=s).digitize(aa) for a,(d,s) in up.items()]
desc={a:d for a,(d,s) in up.items()}
def mafft(names,seqs,out):
    with open("tmp_in.faa","w") as f:
        for n,s in zip(names,seqs): f.write(f">{n}\n{s}\n")
    r=subprocess.run(["mafft","--localpair","--maxiterate","1000","--quiet","tmp_in.faa"],capture_output=True,text=True)
    aln=collections.OrderedDict(); cur=None
    for l in r.stdout.splitlines():
        if l.startswith(">"): cur=l[1:].strip(); aln[cur]=""
        else: aln[cur]+=l.strip()
    open(out,"w").write("".join(f">{n}\n{s}\n" for n,s in aln.items()))
    return aln
def show(aln,ref,lo,hi,width=70):
    # print alignment columns covering ref residues lo..hi (1-based)
    r=aln[ref]; cols=[];k=0
    for i,c in enumerate(r):
        if c!="-":
            k+=1
            if lo<=k<=hi: cols.append(i)
    a,b=cols[0],cols[-1]+1
    for s in range(a,b,width):
        for n,seq in aln.items(): print("   %-14s %s"%(n[:14],seq[s:min(b,s+width)]))
        print()
print("=========== 1. jackhmmer iteration trace, Kp_MATa2")
q=TextSequence(name=b"Kp_MATa2",sequence=S["Kp_MATa2"]).digitize(aa)
markers={"P33417":"Ixr1","Q02486":"Abf2","P25042":"Rox1","Q59PR9":"Hmo1","Q5A220":"Rfg1","Q9HG13":"KlMATA2","Q71U11":"CaMTLA2"}
seenit={}
prev=set()
for it in pyhmmer.jackhmmer(q,pool,max_iterations=6,incE=1e-3,E=1e-3):
    inc={str(h.name) for h in it.hits if h.included}
    new=inc-prev
    print(" iter",it.iteration,"included",len(inc),"new",len(new))
    ex=[(str(h.name),round(h.score,1),desc.get(str(h.name),"")[:60]) for h in it.hits if h.included and str(h.name) in new][:8]
    for e in ex: print("     new:",e)
    for m,n in markers.items():
        if m in inc and m not in seenit: seenit[m]=it.iteration
    prev=inc
print(" first iteration included:",{markers[m]:i for m,i in seenit.items()})
print(" phmmer (iter 1) top hits of Kp_MATa2 with E<=10:")
for hits in pyhmmer.phmmer([q],pool,E=10):
    for h in list(hits)[:12]: print("    %s %.1f %.0e %s"%(h.name,h.score,h.evalue,desc.get(str(h.name),"")[:60]))
print("=========== 2. alignment, MATa2 orthologues and HMG references (region of Kp MATa2 30-130)")
names=["Kp_MATa2","Kl_MATA2","Ca_MTLA2","Dh_MTLA2","Cl_MTLA2","Geo_MATA","Yl_mata2","ctl_NHP6A","ctl_SRY","ctl_IXR1"]
names=[n for n in names if n in S]
pez=[k for k in S if k.startswith("Pez_MAT1-2-1_746128")][0]
names.append(pez)
aln=mafft(names,[S[n] for n in names],"a2.aln.faa")
show(aln,"Kp_MATa2",30,130)
print("=========== 3. ORF checks")
fr=SeqIO.read("FR839631.gb","genbank").seq
# MATa2 minus strand 2740..3276 (1-based) -> start codon at 3276 on minus strand
g=fr[2739:3276].reverse_complement()
print(" MATa2 CDS translates:",str(g.translate())[:30],"... len",len(g)//3-1,"aa")
up_seq=fr[3276:3276+600].reverse_complement()  # upstream region (5') = rc of fr[3276:3876]
ups=str(fr[3276:3876].reverse_complement().translate())
print(" upstream in-frame (5'->3', ends at the ATG):",ups[-60:])
# last stop upstream
i=ups.rfind("*"); print(" last upstream in-frame stop at aa",len(ups)-i,"before ATG; Ms between:",[j for j,c in enumerate(ups) if c=='M' and j>i][:6], "n codons between stop and ATG",len(ups)-i-1)
# downstream: stop codon present
print(" stop codon:",str(g[-3:]))
# MATa1 plus 4058..4735
a1=fr[4057:4735]; print(" MATa1 translate",str(a1.translate())[:20],len(a1)//3-1,"aa; stop",str(a1[-3:]))
u1=str(fr[4057-600:4057].translate()); i=u1.rfind("*"); print(" MATa1 upstream last stop",len(u1)-i,"aa before ATG; Ms in between:",[j-i for j,c in enumerate(u1) if c=='M' and j>i][:6])
print("=========== 4. tblastn of Kp MATa1/MATa2 across Komagataella genomes")
D="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
gen=["GCA_023013455.1_ASM2301345v1","GCA_023013825.1_ASM2301382v1","GCA_023013855.1_ASM2301385v1","GCA_023014715.1_ASM2301471v1","GCA_023013895.1_ASM2301389v1","GCA_023014845.1_ASM2301484v1","GCA_030557195.1_ASM3055719v1","GCA_030567415.1_ASM3056741v1","GCA_000182785.1_ASM18278v1","GCA_900235035.2_KP_7435-4","GCA_000223565.1_PicPas_Mar2011","GCA_001708105.1_ASM170810v1","GCA_023213065.1_ASM2321306v1","GCA_900010155.1_JPS733","GCF_000027005.1_ASM2700v1"]
with open("q_a.faa","w") as f:
    for n in ["Kp_MATa1","Kp_MATa2","Kp_MATalpha1","Kp_MATalpha2"]: f.write(f">{n}\n{S[n]}\n")
for gn in gen:
    tmp=tempfile.mkdtemp(); fa=tmp+"/g.fa"
    with open(fa,"w") as o: o.write(gzip.open(D+gn+".fa.gz","rt").read())
    subprocess.run(["makeblastdb","-in",fa,"-dbtype","nucl","-out",tmp+"/g"],capture_output=True)
    r=subprocess.run(["tblastn","-query","q_a.faa","-db",tmp+"/g","-outfmt","6 qseqid sseqid pident length qstart qend sstart send evalue","-evalue","1e-10","-max_hsps","3"],capture_output=True,text=True).stdout
    best={}
    for l in r.splitlines():
        p=l.split("\t"); q=p[0]
        if q not in best or float(p[3])>float(best[q][3]): best[q]=p
    out=[]
    for q in ["Kp_MATa1","Kp_MATa2","Kp_MATalpha1","Kp_MATalpha2"]:
        if q in best:
            p=best[q]; ql=len(S[q]); out.append("%s %s id%.1f aln%d/%d q%s-%s"%(q[3:],p[1],float(p[2]),int(p[3]),ql,p[4],p[5]))
        else: out.append(q[3:]+" -")
    print(" ",gn[:22]," ;; ".join(out))
    subprocess.run(["rm","-rf",tmp])
