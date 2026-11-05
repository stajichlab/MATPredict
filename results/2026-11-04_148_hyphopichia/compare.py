import glob,yaml,collections,csv
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-04_148_hyphopichia/"
RANK={"low":0,"medium":1,"high":2}
def load(run):
    d={}
    for f in glob.glob(O+run+"/runs/*/detection_report.yaml"):
        r=yaml.safe_load(open(f)); d[f.split("/")[-2]]=[dict(f=l["family"].split(":")[1],idio=l["idiomorph"],conf=l["confidence"],contig=l["contig"],s=l["start"],e=l["end"],genes=sorted({g["gene"] for g in l["gene_evidence"] if not g.get("lost_to")}),
          gev={g["gene"]:(g["identity"],g["coverage"],g["bitscore"],g["status"]) for g in l["gene_evidence"] if not g.get("lost_to")}) for l in r["detected"]]
    return d
b=load("base"); c=load("cand"); print(len(b),len(c))
up=down=same=0; rows=[]
for g in sorted(b):
    x=b[g]; y=c.get(g,[])
    cb=sorted(RANK[l["conf"]] for l in x); cc=sorted(RANK[l["conf"]] for l in y)
    key=lambda L:[(l["f"],l["idio"],l["contig"],l["s"],l["e"],l["conf"],tuple(l["genes"])) for l in L]
    if key(x)==key(y): same+=1; continue
    d="up" if sum(cc)>sum(cb) else "down" if sum(cc)<sum(cb) else "other"
    rows.append((g,d,x,y))
cnt=collections.Counter(r[1] for r in rows); print("same",same,dict(cnt))
for g,d,x,y in rows:
    if d in ("down","other"):
        print(d,g[:16])
        for l in x: print("   base",l["idio"],l["conf"],l["genes"],{k:(round(v[0] or 0),round(v[1] or 0),round(v[2] or 0),v[3][:12]) for k,v in l["gev"].items() if k.startswith("MTLA")})
        for l in y: print("   cand",l["idio"],l["conf"],l["genes"],{k:(round(v[0] or 0),round(v[1] or 0),round(v[2] or 0),v[3][:12]) for k,v in l["gev"].items() if k.startswith("MTLA")})
