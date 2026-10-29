import glob,gzip,os,re,subprocess,sys
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
LIB="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
MP="/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/miniprot"
base={os.path.basename(p):p for p in glob.glob(R+"/2026-10-26_v070_ascomycota/wave_*/small_*/runs/*")}
G=[l.split("\t")[1] for l in open(R+"/2026-10-29_yeast_routing/detail.tsv") if l.startswith("mtl_only")]
WIN=60000
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-29_yeast_routing/flank"
os.makedirs(W,exist_ok=True)
WF=W+"/win.fa"
RF=W+"/flank_refs.faa"
def fasta(path):
    n=None;s=[]
    for l in gzip.open(path,"rt"):
        if l.startswith(">"):
            if n: yield n,"".join(s)
            n=l[1:].split()[0];s=[]
        else: s.append(l.strip())
    if n: yield n,"".join(s)
print("genome\tcontig\tcontig_len\tlocus\tdist_left\tdist_right\tflank_hits")
for g in G:
    for l in open(base[g]+"/detected_loci.gff3"):
        if l.startswith("#"): continue
        c=l.split("\t")
        if c[2]=="MAT_locus": ctg,ls,le=c[0],int(c[3]),int(c[4]);break
    seq=None
    for n,s in fasta(f"{LIB}/{g}.fa.gz"):
        if n==ctg: seq=s;break
    L=len(seq); a=max(0,ls-WIN); b=min(L,le+WIN)
    open(WF,"w").write(">w\n"+seq[a:b]+"\n")
    out=subprocess.run([MP,"--gff","-t2",WF,RF],capture_output=True,text=True); err=out.stderr; out=out.stdout
    if not out: print("#miniprot empty",err[:200],file=sys.stderr)
    hits={}
    for ln in out.splitlines():
        f=ln.split("\t")
        if len(f)>8 and f[2]=="mRNA":
            at=dict(x.split("=",1) for x in f[8].split(";") if "=" in x)
            nm=re.search(r"name=(\w+)",ln); 
            nm=at.get("Target","").split()[0]
            m=re.search(r"\|name=(\w+)",at.get("Target","")+" "+ln)
            gname=("PAP1" if "PAP1" in ln else "OBP1" if "OBP1" in ln else "PIK1" if "PIK1" in ln else "?")
            ident=float(at.get("Identity",0))*100
            rec=(ident,int(f[3])+a,int(f[4])+a,f[6])
            if gname not in hits or rec[0]>hits[gname][0]: hits[gname]=rec
    hs=";".join(f"{k}:{v[0]:.0f}%@{v[1]}-{v[2]}{v[3]}" for k,v in sorted(hits.items()))
    print(f"{g}\t{ctg}\t{L}\t{ls}-{le}\t{ls-1}\t{L-le}\t{hs}")
