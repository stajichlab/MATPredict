import gzip,subprocess,os,sys,tempfile,shutil
from multiprocessing import Pool
B="/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/"
G="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
W=os.path.abspath("chrIII_windows.fa")
TMP=os.path.abspath("tmpdb")
def run(a):
    d=os.path.join(TMP,a); os.makedirs(d,exist_ok=True)
    try:
        seq={};k=None;buf=[]
        with gzip.open(G+a+".fa.gz","rt") as f, open(d+"/g.fa","w") as o:
            for l in f:
                o.write(l)
                if l[0]=='>':
                    if k: seq[k]="".join(buf)
                    k=l[1:].split()[0];buf=[]
                else: buf.append(l.strip())
            if k: seq[k]="".join(buf)
        subprocess.run([B+"makeblastdb","-in",d+"/g.fa","-dbtype","nucl"],capture_output=True,check=True)
        r=subprocess.run([B+"blastn","-query",W,"-db",d+"/g.fa","-outfmt","6 qseqid sseqid pident length qstart qend sstart send","-evalue","1e-50","-max_hsps","5","-max_target_seqs","10"],capture_output=True,text=True,check=True).stdout
        hits={}
        for l in r.splitlines():
            q,s,pid,ln,qs,qe,ss,se=l.split("\t"); ln=int(ln);qs,qe,ss,se=map(int,(qs,qe,ss,se))
            if float(pid)<90: continue
            hits.setdefault(q.split("_")[0],[]).append((ln,s,qs,qe,ss,se))
        out=[a]
        for c in ["HML","MAT","HMR"]:
            L=hits.get(c+"left",[]); R=hits.get(c+"right",[])
            # left flank: hit covering the flank's inner end (qend >= 8500 of 9000 / 8001 windows) ; take the hit with max qend then longest
            Lin=[h for h in L if h[3]>=len_w[c+"left"]-300]
            Rin=[h for h in R if h[2]<=300]
            if not Lin or not Rin:
                out.append(f"{c}:flank_missing(L={len(Lin)},R={len(Rin)})"); continue
            lh=max(Lin,key=lambda h:h[0]); rh=max(Rin,key=lambda h:h[0])
            if lh[1]!=rh[1]:
                # maybe several; try any pair same contig
                pairs=[(x,y) for x in Lin for y in Rin if x[1]==y[1]]
                if not pairs:
                    # edge distance of left inner end on its contig
                    out.append(f"{c}:split_contigs"); continue
                lh,rh=max(pairs,key=lambda p:p[0][0]+p[1][0])
            lp=lh[5]; rp=rh[4]  # genomic coords of left-inner-end and right-inner-start
            lo,hi=sorted((lp,rp))
            if hi-lo>50000: out.append(f"{c}:far({hi-lo})"); continue
            s=seq[lh[1]][lo:hi-1]
            out.append(f"{c}:gap={hi-lo-1};N={s.upper().count('N')};contig={lh[1]};pos={lo}")
        return "\t".join(out)
    except Exception as e:
        return f"{a}\tERROR {e}"
    finally:
        shutil.rmtree(d,ignore_errors=True)
len_w={}
k=None
for l in open(W):
    if l[0]=='>': k=l[1:].strip().split("_")[0]
    else: len_w[k]=len(l.strip())
if __name__=="__main__":
    ids=[l.strip() for l in open(sys.argv[1])]
    with Pool(4) as p, open(sys.argv[2],"w") as o:
        for r in p.imap_unordered(run,ids): o.write(r+"\n"); o.flush()
