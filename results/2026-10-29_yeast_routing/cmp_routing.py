import json,glob,os,sys,collections
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
W=R+"/2026-10-29_yeast_routing"
groups=json.load(open(W+"/yeast_groups.json"))
base={}
for p in glob.glob(R+"/2026-10-26_v070_ascomycota/wave_*/small_*/runs/*"):
    base[os.path.basename(p)]=p
def loci(d):
    f=d+"/detected_loci.gff3"
    out=[]
    if not os.path.exists(f): return None
    for l in open(f):
        if l.startswith("#"): continue
        c=l.rstrip("\n").split("\t")
        if len(c)>8 and c[2]=="MAT_locus":
            a=dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
            out.append((c[0],int(c[3]),int(c[4]),a["family"].split(":")[-1],a.get("confidence","")))
    return sorted(out)
def fams(L): return sorted(x[3] for x in L) if L is not None else None
res=collections.defaultdict(lambda: collections.Counter())
detail=[]
for grp,gl in groups.items():
    for g in gl:
        b=loci(base[g]) if g in base else None
        v1=loci(W+"/out_v1/runs/"+g); v2=loci(W+"/out_v2/runs/"+g)
        fb,f1,f2=fams(b),fams(v1),fams(v2)
        res[grp]["n"]+=1
        res[grp]["v1_same_as_base"]+= f1==fb
        res[grp]["v2_same_as_base"]+= f2==fb
        res[grp]["v1_eq_v2"]+= f1==f2
        res[grp]["base_dupe"]+= b is not None and len(b)>1
        res[grp]["v1_dupe"]+= v1 is not None and len(v1)>1
        res[grp]["v2_dupe"]+= v2 is not None and len(v2)>1
        res[grp]["base_none"]+= not b
        res[grp]["v1_none"]+= not v1
        res[grp]["v2_none"]+= not v2
        detail.append((grp,g,fb,f1,f2))
for g,c in res.items(): print(g,dict(c))
with open(sys.argv[1],"w") as o:
    for d in detail: o.write("\t".join(str(x) for x in d)+"\n")
