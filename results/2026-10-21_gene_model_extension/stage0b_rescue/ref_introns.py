"""Intron positions (aa, along the protein) of every curated reference. python ref_introns.py DB_ROOT OUT_TSV"""
import glob, yaml, sys, csv
rows=[]
for md in sorted(glob.glob(f"{sys.argv[1]}/*/*/*/metadata.yaml")):
    if "/candidates/" in md or "/classifiers/" in md: continue
    d=yaml.safe_load(open(md))
    for g in d.get("genes") or []:
        if g["name"] not in ("MAT1-1-1","MAT1-2-1") or g.get("role")!="core_MAT" or not g.get("exons"): continue
        ex=[(e["start"],e["end"]) for e in g["exons"]]
        # exons are listed in transcript order; CDS intron boundaries = cumulative lengths
        cum=0; pos=[]
        for a,b in ex[:-1]:
            cum+=b-a+1; pos.append(round(cum/3,1))
        rows.append(dict(record=d["record_id"],gene=g["name"],strand=g["strand"],nexons=len(ex),cds_aa=round(sum(b-a+1 for a,b in ex)/3,1),intron_aa=";".join(map(str,pos))))
with open(sys.argv[2],"w",newline="") as f:
    w=csv.DictWriter(f,rows[0].keys(),delimiter="\t"); w.writeheader(); w.writerows(rows)
print(len(rows),"references with exon structure")
