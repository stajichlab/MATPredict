import sys,os,re,csv,subprocess,tempfile,collections,shutil
sys.path.insert(0,"/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit")
import taph_mat_candidates as t
P="/bigdata/stajichlab/jstajich/projects/MATPredict"; BIN=P+"/.pixi/envs/default/bin"; I=P+"/results/sra_reads_typing/wave_inputs/genomes"
O=P+"/results/2026-10-28_taphrinomycotina_mat_candidates"; OUT=O+"/probe"; os.makedirs(OUT,exist_ok=True)
GENE={"matPi":"Pi","mat1-Pi":"Pi","matMc":"Mc","mat1-Mc":"Mc","matMi":"Mi","mat1-Mi":"Mi","matPc":"Pc","mat1-Pc":"Pc"}
refs={}      # group -> [(id,seq)]
for h,s in t.read_fasta(O+"/references.faa"):
    if h.endswith("|PM"):
        g=GENE.get(h.split("|")[1]); refs.setdefault(g,[]).append((h,s))
# candidate pool (strict run): protein sequences with species and best gene
seqs=dict(t.read_fasta(O+"/candidates.faa")); seqs={k.split()[0]:v for k,v in seqs.items()}
cand=[]
for r in csv.DictReader(open(O+"/candidates.tsv"),delimiter="\t"):
    if r["ref_class"]=="PM" and not r["species"].startswith("Schizosaccharomyces pombe") and float(r["identity"])>=0.35 and r["candidate"] in seqs:
        cand.append((r["species"],GENE.get(r["ref_gene"]),float(r["identity"]),seqs[r["candidate"]].rstrip("*"),r["candidate"]))
# windows (genome, contig, start, end, label, genus)
W=[("GCA_001929475.1_NeoIirr1.0","LXFE01000411.1",5949,25065,"Neolecta matMc+matPi window","Neolecta"),
   ("GCA_001929475.1_NeoIirr1.0","LXFE01003167.1",5,851,"Neolecta matMc contig-edge","Neolecta"),
   ("GCA_001929475.1_NeoIirr1.0","LXFE01003595.1",11513,12163,"Neolecta matPi","Neolecta"),
   ("GCA_051599425.1_Sai_col_YB2330_v1","AP039747.1",177119,198616,"S. coloradoensis multi-gene window","Saitoella"),
   ("GCA_051599425.1_Sai_col_YB2330_v1","AP039745.1",105707,106351,"S. coloradoensis matMc-like","Saitoella"),
   ("GCF_001661265.1_Saico1","NW_017567017.1",185699,210006,"S. complicata multi-gene window","Saitoella")]
# controls: detected candidates of Protomyces / Taphrina species (moderate identity)
for r in csv.DictReader(open(O+"/candidates.tsv"),delimiter="\t"):
    if r["species"].split()[0] in ("Protomyces","Taphrina") and r["ref_class"]=="PM" and "detected" in r["report_locus"] and len([w for w in W if w[5]==r["species"].split()[0]])<3:
        W.append((r["genome"],r["contig"],int(r["start"]),int(r["end"]),f"control {r['species']} {r['ref_gene']} id{r['identity']}",r["species"].split()[0]))
def mafft(seqlist):
    d=tempfile.mkdtemp(dir="/tmp"); open(d+"/a.faa","w").write("".join(f">{i}\n{s}\n" for i,s in seqlist))
    p=subprocess.run([BIN+"/mafft","--auto","--quiet",d+"/a.faa"],capture_output=True,text=True)
    aln={};n=None
    for l in p.stdout.splitlines():
        if l.startswith(">"): n=l[1:]; aln[n]=[]
        elif n: aln[n].append(l.strip())
    return {k:"".join(v) for k,v in aln.items()}
def consensus(aln,minfrac=0.4):
    L=len(next(iter(aln.values()))); out=[]
    for i in range(L):
        col=[a[i] for a in aln.values() if a[i]!="-"]
        if len(col)>=minfrac*len(aln): out.append(collections.Counter(col).most_common(1)[0][0])
    return "".join(out)
def hmm_from(aln,name):
    import pyhmmer
    alpha=pyhmmer.easel.Alphabet.amino()
    msa=pyhmmer.easel.TextMSA(name=name.encode(),sequences=[pyhmmer.easel.TextSequence(name=k.encode(),sequence=v) for k,v in aln.items()]).digitize(alpha)
    b=pyhmmer.plan7.Builder(alpha); bg=pyhmmer.plan7.Background(alpha)
    hmm,_,_=b.build_msa(msa,bg); return hmm,alpha,bg
def sixframe(dna):
    from Bio.Seq import Seq
    out=[]
    for strand,sq in (("+",Seq(dna)),("-",Seq(dna).reverse_complement())):
        for f in range(3):
            sub=str(sq[f:]); sub=sub[:len(sub)//3*3]
            out.append((strand,f,str(Seq(sub).translate())))
    return out
def hmm_scan(hmm,alpha,bg,dna):
    import pyhmmer
    best=None
    for strand,f,prot in sixframe(dna):
        for chunk in prot.split("*"):
            if len(chunk)<30: continue
            sq=pyhmmer.easel.TextSequence(name=b"x",sequence=chunk).digitize(alpha)
            for hits in pyhmmer.hmmsearch([hmm],[sq],background=bg,E=10):
                for h in hits:
                    if best is None or h.evalue<best[0]: best=(h.evalue,h.score,len(chunk))
    return best
def run_tools(window_fa,query_fa,maxintron=5000):
    res={}
    p=subprocess.run([BIN+"/miniprot","-G",str(maxintron),"--gff",window_fa,query_fa],capture_output=True,text=True)
    mp=t.parse_gff(p.stdout); ql={h:len(s) for h,s in t.read_fasta(query_fa)}
    mp=[(m["ident"],(m["tend"]-m["tstart"]+1)/ql.get(m["target"],1),len(m["cds"]),m["score"]) for m in mp]
    res["miniprot"]=max(mp,key=lambda x:x[0]*x[1]) if mp else None
    p=subprocess.run([BIN+"/exonerate","--model","protein2genome","--query",query_fa,"--target",window_fa,"--maxintron",str(maxintron),"--refine","region","--showtargetgff","yes","--showalignment","no"],capture_output=True,text=True)
    al=[];cur=None
    for ln in p.stdout.splitlines():
        f=ln.split("\t")
        if len(f)<9 or ln.startswith("#"): continue
        if f[2]=="gene":
            m=re.search(r"identity ([\d.]+)",f[8]); q=re.search(r"sequence (\S+)",f[8]); cur=dict(score=float(f[5]),ident=float(m.group(1))/100 if m else 0,q=q.group(1) if q else "",cds=0,exons=0); al.append(cur)
        elif f[2]=="cds" and cur is not None: cur["cds"]+=int(f[4])-int(f[3])+1
        elif f[2]=="exon" and cur is not None: cur["exons"]+=1
    ex=[(a["ident"],a["cds"]/3/ql.get(a["q"],1) if a["q"] in ql else 0,a["exons"],a["score"]) for a in al]
    res["exonerate"]=max(ex,key=lambda x:x[0]*x[1]) if ex else None
    return res
fmt=lambda r:"-" if not r else "id%.2f cov%.2f ex%d s%.0f"%r
print("windows:",len(W)); rows=[]
for g,c,s,e,label,genus in W:
    seqs_g=t.load_genome(f"{I}/{g}.fa.gz"); s0=max(1,s-1500); win=seqs_g[c][s0-1:e+1500]
    d=tempfile.mkdtemp(dir="/tmp"); open(d+"/w.fa","w").write(f">{c}\n{win}\n")
    print("\n##",label,f"{c}:{s}-{e}")
    for grp in ("Pi","Mc","Mi"):
        base=refs.get(grp,[])
        extra=[(f"{x[0].replace(chr(32),chr(95))}|{x[4]}",x[3]) for x in cand if x[1]==grp and x[0].split()[0]!=genus]
        # dedupe
        seen=set(); extra=[(i,q) for i,q in extra if not (q in seen or seen.add(q))]
        sets={"A refs (pipeline)":base,"B refs+close reps":base+extra}
        aln=mafft(base+extra) if len(base+extra)>=3 else None
        cons=consensus(aln) if aln else None
        if cons and len(cons)>=40: sets["C consensus"]=[(f"cons_{grp}",cons)]
        for name,ql_ in sets.items():
            qf=d+"/q.faa"; open(qf,"w").write("".join(f">{i}\n{q}\n" for i,q in ql_))
            r=run_tools(d+"/w.fa",qf)
            print(f"  {grp:2s} {name:19s} n_query={len(ql_):2d}  miniprot {fmt(r['miniprot']):28s} exonerate {fmt(r['exonerate'])}")
            rows.append((label,grp,name,len(ql_),fmt(r["miniprot"]),fmt(r["exonerate"])))
        if aln:
            hmm,alpha,bg=hmm_from(aln,grp)
            b=hmm_scan(hmm,alpha,bg,win)
            print(f"  {grp:2s} HMM ({len(aln)} seqs, M={hmm.M}): best 6-frame hit "+("-" if not b else "E=%.1e score %.1f in %d aa"%b))
            rows.append((label,grp,"HMM",len(aln),"-" if not b else "E=%.1e bits=%.1f"%(b[0],b[1]),""))
with open(OUT+"/probe_table.tsv","w") as fh:
    w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["window","gene","query_set","n_queries","miniprot","exonerate"]); w.writerows(rows)
