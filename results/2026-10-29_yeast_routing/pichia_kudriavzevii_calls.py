import glob,os
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
for g in ["GCF_003054445.1_ASM305444v1","GCA_003054405.1_ASM305440v1","GCA_003033855.1_ASM303385v1","GCA_014873065.1_ASM1487306v1"]:
    p=glob.glob(R+"/2026-10-26_v070_ascomycota/wave_*/small_*/runs/"+g)[0]
    print("##",g)
    for l in open(p+"/detected_loci.gff3"):
        if l.startswith("#"): continue
        c=l.rstrip("\n").split("\t"); a=dict(x.split("=",1) for x in c[8].split(";") if "=" in x)
        if c[2]=="MAT_locus": print(" LOCUS",c[0],c[3],c[4],a["family"],a["confidence"],a.get("idiomorph"),"frag="+a["fragmented"],a.get("reference_records","")[:80])
        elif c[2]=="gene" and a.get("present")=="true": print("    ",a["Name"],a.get("role"),c[3],c[4],c[6],a.get("identity","")[:5],a.get("status",""))
