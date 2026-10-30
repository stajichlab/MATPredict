import glob,os,collections,json
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-29_yeast_routing/"
rows=[l.rstrip("\n").split("\t") for l in open(R+"detail.tsv")]
G=sorted({r[1] for r in rows if r[0]=="ogataea_other" or "MATog" in r[4]})
out=[]
for g in G:
    loci=[];cur=None
    for l in open(R+"out_v2/runs/"+g+"/detected_loci.gff3"):
        if l.startswith("#"): continue
        c=l.rstrip("\n").split("\t"); a=dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
        if c[2]=="MAT_locus": cur={"contig":c[0],"start":int(c[3]),"end":int(c[4]),"fam":a["family"].split(":")[1],"idio":a.get("idiomorph"),"conf":a["confidence"],"genes":[]}; loci.append(cur)
        elif c[2]=="gene" and a.get("present")=="true" and cur: cur["genes"].append(a["Name"])
    loci=[l for l in loci if l["fam"]=="MATog"]
    if not loci: pattern="no_MATog"
    elif len(loci)==1: pattern="one:"+loci[0]["idio"]+":"+"+".join(sorted(set(loci[0]["genes"])))
    else:
        ids=sorted(l["idio"] for l in loci)
        same="same_contig" if len({l["contig"] for l in loci})==1 else "diff_contigs"
        pattern=f"{len(loci)}:{'+'.join(ids)}:{same}"
    out.append((g,pattern))
c=collections.Counter(p.split(":",1)[0] if p.startswith("one") else p for g,p in out)
for k,v in c.most_common(): print(v,k)
print("--one-locus gene sets")
print(collections.Counter(p for g,p in out if p.startswith("one")).most_common(10))
print("--edge proxy (locus starts within 500 bp of its contig start) among the 71 a+alpha genomes")
n=0;e=collections.Counter()
for g,p in out:
    if p!="2:MATa+MATalpha:diff_contigs": continue
    loci=[];
    for l in open(R+"out_v2/runs/"+g+"/detected_loci.gff3"):
        if l.startswith("#"): continue
        c=l.split("\t")
        if c[2]=="MAT_locus" and "MATog" in c[8]: loci.append(int(c[3]))
    n+=1; e[sum(1 for s in loci if s<=500)]+=1
print(n,dict(e))
