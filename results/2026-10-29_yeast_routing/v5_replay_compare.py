import glob,os,json,collections,re
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
W=R+"/2026-10-29_yeast_routing"
groups=json.load(open(W+"/yeast_groups.json"))
base={os.path.basename(p):p for p in glob.glob(R+"/2026-10-26_v070_ascomycota/wave_*/small_*/runs/*")}
def loci(d):
    f=d+"/detected_loci.gff3"
    if not os.path.exists(f): return None
    out=[]
    for l in open(f):
        if l.startswith("#"): continue
        c=l.rstrip("\n").split("\t")
        if len(c)>8 and c[2]=="MAT_locus":
            a=dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
            out.append((c[0],int(c[3]),int(c[4]),a["family"].split(":")[-1],a["confidence"],a.get("family_status",""),a.get("cross_call","")))
    return out
res=collections.defaultdict(collections.Counter); st=collections.Counter(); lost=[]
for grp,gl in groups.items():
    for g in gl:
        b=loci(base[g]); v=loci(W+"/out_v5/runs/"+g)
        if v is None: res[grp]["missing"]+=1; continue
        res[grp]["n"]+=1
        res[grp]["base_2plus"]+=len(b)>1; res[grp]["v3_2plus"]+=len(v)>1
        res[grp]["base_none"]+=not b; res[grp]["v3_none"]+=not v
        res[grp]["conf_changed"]+= [x[4] for x in b]!=[x[4] for x in v] and len(b)==len(v)
        for x in v:
            if x[5]: st[(grp,x[5])]+=1
        if grp=="mtl_only" and not v: lost.append(g)
for g,c in res.items(): print(g,dict(c))
print(dict(st)); print("mtl_only lost:",lost)
print("--cross181 v3 patterns")
pat=collections.Counter(); ov=[]
for g in groups["cross181"]:
    v=loci(W+"/out_v5/runs/"+g)
    pat[tuple(sorted(x[3] for x in v))]+=1
    fams={x[3] for x in v}
    if len(fams)>1:
        a=[x for x in v if x[3]=="MATsc"][0]; m=[x for x in v if x[3]=="MTL"][0]
        sh=min(a[2],m[2])-max(a[1],m[1])+1 if a[0]==m[0] else -1
        ov.append((a[0]==m[0], sh))
for k,n in pat.most_common(): print(n,k)
print("MATsc+MTL survivors:",len(ov),"same contig:",sum(1 for o in ov if o[0]),"overlap>0:",sum(1 for o in ov if o[1]>0))
