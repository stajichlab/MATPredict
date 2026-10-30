import glob,os,json,collections
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"; W=R+"/2026-10-29_yeast_routing"
groups=json.load(open(W+"/yeast_groups.json"))
tx={l.split("\t")[0]:l.split("\t")[1].strip() for l in open(W+"/yeast_list.tsv")}
base={os.path.basename(p):p for p in glob.glob(R+"/2026-10-26_v070_ascomycota/wave_*/small_*/runs/*")}
def loci(d):
    out=[]
    for l in open(d+"/detected_loci.gff3"):
        if l.startswith("#"): continue
        c=l.rstrip("\n").split("\t")
        if c[2]=="MAT_locus":
            a=dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
            out.append((a["family"].split(":")[1],a["confidence"],a.get("family_status",""),a.get("cassette","")))
    return out
# the 25 MTL-wins genomes: which family survives now
tie=[g for g in groups["cross181"] if {x[0] for x in loci(W+"/out_v3/runs/"+g)}=={"MTL"}]
c=collections.Counter()
for g in tie:
    v=loci(W+"/out_v5/runs/"+g); c[tuple(sorted({(x[0],x[2]) for x in v}))]+=1
print("25 formerly MTL-wins:",len(tie)); [print(" ",n,k) for k,n in c.most_common()]
# status counts over cross181 in v5
s=collections.Counter()
for g in groups["cross181"]:
    for x in loci(W+"/out_v5/runs/"+g): s[(x[0],x[2])]+=1
print("cross181 v5 loci by family/status:",dict(s))
# lost call detail
g="GCA_019232045.1_Wolfe_Odeg_UCD465"; print(g,"v5:",loci(W+"/out_v5/runs/"+g),"v070:",loci(base[g]))
# conf changes: which are family-changed Ogataea
oc=0;other=[]
for grp in ["cross181","sc_only40","ogataea_other"]:
    for g in groups[grp]:
        b=loci(base[g]); v=loci(W+"/out_v5/runs/"+g)
        if len(b)==len(v) and [x[1] for x in b]!=[x[1] for x in v]:
            if any(x[0]=="MATog" for x in v): oc+=1
            else: other.append(g)
print("conf changed, Ogataea (MATog now):",oc,"other:",other)
