import sys
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-29_yeast_routing/out_v2/runs/"
for g in sys.argv[1:]:
    print("##",g)
    cur=None
    for l in open(R+g+"/detected_loci.gff3"):
        if l.startswith("#"): continue
        c=l.rstrip("\n").split("\t"); a=dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
        if c[2]=="MAT_locus": print(" LOCUS",c[0],c[3],c[4],a["family"],a["confidence"],a.get("idiomorph"),"frag="+a["fragmented"])
        elif c[2]=="gene" and a.get("present")=="true": print("   ",a["Name"],a.get("role"),c[3],c[4],c[6],a.get("identity","")[:5],a.get("status",""))
