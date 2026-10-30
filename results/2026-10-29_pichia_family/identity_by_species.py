
import sys,collections,os
R=sys.argv[1]+"/out_pichia_after/runs/"
names=dict(l.rstrip("\n").split("\t") for l in open(sys.argv[1]+"/pichia_names.tsv"))
agg=collections.defaultdict(list)
for g in sorted(os.listdir(R)):
    f=R+g+"/detected_loci.gff3"
    if not os.path.exists(f): continue
    sp=names.get(g,"?")
    cur=None
    for l in open(f):
        if l.startswith("#"): continue
        c=l.rstrip("\n").split("\t"); a=dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
        if c[2]=="MAT_locus": cur=[sp,g,a.get("idiomorph"),a.get("confidence"),a.get("locus_class"),[]]; agg[sp].append(cur)
        elif c[2]=="gene" and a.get("present")=="true" and a.get("role")=="core_MAT" and cur: cur[5].append((a["Name"],a.get("identity","")[:5]))
for sp,ls in sorted(agg.items(), key=lambda kv:-len(kv[1])):
    core=[x for x in ls if x[5]]
    ids=[float(i) for x in ls for n,i in x[5] if i]
    print(sp, "loci",len(ls),"with core",len(core), "core identity min/median/max", (min(ids),sorted(ids)[len(ids)//2],max(ids)) if ids else None, collections.Counter(x[3] for x in ls))

