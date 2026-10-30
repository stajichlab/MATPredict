import glob,os,json
try:
    import yaml
except Exception: yaml=None
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"; W=R+"/2026-10-29_yeast_routing"
groups=json.load(open(W+"/yeast_groups.json"))
base={os.path.basename(p):p for p in glob.glob(R+"/2026-10-26_v070_ascomycota/wave_*/small_*/runs/*")}
def gff(d):
    L=[];cur=None
    for l in open(d+"/detected_loci.gff3"):
        if l.startswith("#"): continue
        c=l.rstrip("\n").split("\t"); a=dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
        if c[2]=="MAT_locus": cur={"contig":c[0],"s":int(c[3]),"e":int(c[4]),"fam":a["family"].split(":")[1],"conf":a["confidence"],"idio":a.get("idiomorph"),"genes":[]};L.append(cur)
        elif c[2]=="gene" and a.get("present")=="true" and cur: cur["genes"].append(f'{a["Name"]}:{a.get("identity","")[:4]}:{a.get("status","").replace("polished_","")}')
    return L
n=0
for g in groups["cross181"]:
    v3=gff(W+"/out_v3/runs/"+g)
    if {x["fam"] for x in v3}!={"MTL"}: continue
    n+=1; b=gff(base[g])
    print("##",g)
    for x in b:
        print(f'  {x["fam"]:6} {x["conf"]:6} {x["idio"]:8} {x["contig"]}:{x["s"]}-{x["e"]}  '+" ".join(x["genes"]))
print(n)
