import csv, glob, os, collections, yaml, sys, re
L=getattr(yaml,"CSafeLoader",yaml.SafeLoader)
sc={(r["genome"],r["contig"],r["start"]):r for r in csv.DictReader(open(sys.argv[2]),delimiter="\t")}
def idm(n):
    m=re.match(r"MAT1-([12])-",n); return None if not m else "MAT1-"+m.group(1)
def vote(ev, drop):
    best={}
    for e in ev:
        n=e["gene"]; i=idm(n)
        if not i or n in drop or e.get("lost_to"): continue
        s=e.get("bitscore") if e.get("bitscore") is not None else e.get("identity")
        if s is not None and s>best.get(i,-1): best[i]=s
    if not best: return None
    r=sorted(best.items(), key=lambda kv:-kv[1])
    if len(r)==1 or r[0][1]>r[1][1]: return r[0][0]
    return "undetermined"
variants={"all (replays campaign)":set(), "drop MAT1-1-2":{"MAT1-1-2"}, "drop MAT1-1-2 and -3":{"MAT1-1-2","MAT1-1-3"},
          "core only (1-1-1, 1-2-1)":{"MAT1-1-2","MAT1-1-3","MAT1-1-4","MAT1-1-5","MAT1-2-2","MAT1-2-3","MAT1-2-4","MAT1-2-10"}}
tot=0; same_as_campaign=collections.Counter(); change=collections.defaultdict(collections.Counter); prof=collections.defaultdict(collections.Counter)
for p in glob.glob(sys.argv[1]+"/**/detection_report.yaml", recursive=True):
    g=os.path.basename(os.path.dirname(p)); doc=yaml.load(open(p),Loader=L)
    for c in doc["detected"]:
        if c["family"]!="Ascomycota:MAT": continue
        tot+=1
        k=(g,c["contig"],str(c["start"])); s=sc.get(k)
        for name,drop in variants.items():
            v=vote(c["gene_evidence"],drop)
            if name.startswith("all"): same_as_campaign[v==c["idiomorph"]]+=1
            elif v is not None and v!=c["idiomorph"]:
                change[name][(c["idiomorph"],v)]+=1
            if s and v is not None and not name.startswith("all"):
                pc=s["profile_call"]
                prof[name][("old label" if False else "new label agrees profile", v==pc, "old agrees", c["idiomorph"]==pc)]+=1
print("loci",tot,"vote replay equals campaign label:",dict(same_as_campaign))
for n in variants:
    if n.startswith("all"): continue
    print("\n",n,"- loci whose label changes:",sum(change[n].values()),dict(change[n]))
    pp=prof[n]; 
    gained=sum(v for k,v in pp.items() if k[1] and not k[3]); lost=sum(v for k,v in pp.items() if not k[1] and k[3])
    print("   among profile-scored loci (n=%d): agree with profile old %d -> new %d  (gained %d, lost %d)"%(
        sum(pp.values()), sum(v for k,v in pp.items() if k[3]), sum(v for k,v in pp.items() if k[1]), gained, lost))
