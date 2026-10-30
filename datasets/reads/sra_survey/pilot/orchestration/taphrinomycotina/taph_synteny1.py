import re,subprocess,os,collections,sys,csv
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates/synteny"; A=O+"/annot"; DIAMOND="/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/diamond"
def fasta(p):
    d={};n=None;desc={}
    for l in open(p):
        if l.startswith(">"): n=l[1:].split()[0]; d[n]=[]; desc[n]=l[1:].strip()
        elif n: d[n].append(l.strip())
    return {k:"".join(v) for k,v in d.items()},desc
def gff_genes(gff):
    genes={};rna={};cds=collections.defaultdict(list)
    for l in open(gff):
        if l.startswith("#"): continue
        f=l.rstrip("\n").split("\t")
        if len(f)<9: continue
        at=dict(kv.split("=",1) for kv in f[8].split(";") if "=" in kv)
        if f[2]=="gene": genes[at["ID"]]=dict(contig=f[0],start=int(f[3]),end=int(f[4]),strand=f[6],name=at.get("Name",""),id=at["ID"],pid=None)
        elif f[2] in ("mRNA","transcript") and "Parent" in at: rna[at["ID"]]=at["Parent"]
        elif f[2]=="CDS" and "Parent" in at:
            pid=at.get("protein_id") or at.get("ID","").replace("cds-","")
            cds[at["Parent"]].append(pid)
    for r,gid in rna.items():
        if gid in genes and cds.get(r): genes[gid]["pid"]=cds[r][0]
    # gene-level CDS (no mRNA, e.g. JGI-style)
    for gid,g in genes.items():
        if g["pid"] is None and cds.get(gid): g["pid"]=cds[gid][0]
    return genes
def region_genes(label,contig,s,e,pad):
    G=gff_genes(f"{A}/{label}.genomic.gff"); prot,desc=fasta(f"{A}/{label}.protein.faa")
    sel=[g for g in G.values() if g["contig"]==contig and g["end"]>=s-pad and g["start"]<=e+pad and g["pid"] in prot]
    sel.sort(key=lambda g:g["start"]); return sel,prot,desc,G
os.makedirs(O+"/dm",exist_ok=True)
# reference flank set: Pjir +/- 50 kb of the PM locus
sel,prot,desc,G=region_genes("Pjir_RU7","NW_017264783.1",211240,221394,50000)
pm=[g for g in sel if g["end"]>=211240 and g["start"]<=221394]
with open(O+"/pjir_flank.faa","w") as fh:
    for g in sel: fh.write(f">{g['id']}|{g['pid']}|{g['start']}-{g['end']}|{'PM' if g in pm else ('L' if g['end']<211240 else 'R')}\n{prot[g['pid']]}\n")
print("Pjir flank genes:",len(sel),"| PM-locus genes:",[(g['id'],g['start'],g['end'],desc[g['pid']][:70]) for g in pm])
targets={"Pcar_B80":"Pcar_B80","Spombe":"Spombe","Tdef_PYCC":"Tdef_PYCC","Scomp":"Scomp","Scoloradoensis":"Scoloradoensis","Neolecta":"Neolecta","Plactucae":"Plactucae"}
rows=[]
for lab in targets:
    pf=f"{A}/{lab}.protein.faa"
    if not os.path.exists(pf): continue
    db=f"{O}/dm/{lab}"
    subprocess.run([DIAMOND,"makedb","--in",pf,"-d",db,"--quiet"],check=True)
    out=f"{O}/dm/pjir_vs_{lab}.tsv"
    subprocess.run([DIAMOND,"blastp","-q",O+"/pjir_flank.faa","-d",db,"-o",out,"--more-sensitive","-e","1e-5","-k","3","--outfmt","6","qseqid","sseqid","pident","length","evalue","bitscore","--quiet"],check=True)
    Gt=gff_genes(f"{A}/{lab}.genomic.gff"); byp={g["pid"]:g for g in Gt.values() if g["pid"]}
    best={}
    for l in open(out):
        q,s_,pid,ln,ev,bs=l.split("\t")
        if q not in best: best[q]=(s_,float(pid),float(ev),float(bs.strip()))
    for g in sel:
        q=f"{g['id']}|{g['pid']}|{g['start']}-{g['end']}|{'PM' if g in pm else ('L' if g['end']<211240 else 'R')}"
        if q in best:
            sid,pid,ev,bs=best[q]; tg=byp.get(sid)
            if tg: rows.append((lab,g["id"],g["start"],q.split("|")[3],tg["contig"],tg["start"],tg["end"],tg["strand"],pid,ev,sid))
with open(O+"/pjir_flank_orthologs.tsv","w") as fh:
    w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["target","pjir_gene","pjir_start","side","contig","start","end","strand","pident","evalue","protein"]); w.writerows(rows)
# synteny blocks: per target and contig, hits within 150 kb windows
for lab in targets:
    r=[x for x in rows if x[0]==lab]
    by=collections.defaultdict(list)
    for x in r: by[x[4]].append(x)
    print(f"\n== {lab}: {len(r)} of {len(sel)} flank genes have a best hit (e<1e-5)")
    for c,v in sorted(by.items(),key=lambda kv:-len(kv[1]))[:3]:
        v.sort(key=lambda x:x[5]); span=v[-1][6]-v[0][5]
        print(f"  contig {c}: {len(v)} hits spanning {span//1000} kb ({v[0][5]}-{v[-1][6]}); sides {''.join(x[3][0] for x in v)}  pjir order(start): {[x[2]//1000 for x in v][:20]}")
