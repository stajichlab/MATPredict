import glob,os,collections
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-29_yeast_routing/out_v4/runs/"
c=collections.Counter();ex=[]
for d in sorted(glob.glob(R+"*")):
    L=[]
    for l in open(d+"/detected_loci.gff3"):
        if l.startswith("#"): continue
        f=l.rstrip("\n").split("\t")
        if f[2]=="MAT_locus":
            a=dict(x.split("=",1) for x in f[8].split(";") if "=" in x); L.append((a["family"].split(":")[1],a.get("idiomorph"),a.get("cassette",""),a.get("family_status","")))
    og=[x for x in L if x[0]=="MATog"]
    key=(len(og),tuple(sorted(x[2] for x in og)))
    c[key]+=1
    if len(ex)<3 and len(og)==2: ex.append((os.path.basename(d),L))
for k,v in c.most_common(): print(v,k)
for e in ex: print(e)
print("other-family loci:",sum(1 for d in glob.glob(R+"*") for l in open(d+"/detected_loci.gff3") if "\tMAT_locus\t" in l and "MATog" not in l))
