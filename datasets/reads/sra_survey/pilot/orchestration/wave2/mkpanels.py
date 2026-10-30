import gzip, sys, os
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-b317ece/src")
from MATPredict.reads.panel import Panel
D="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
OUT="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/panels/"
SPEC={
 "pyricularia_oryzae":{"MAT1-1":("GCA_043231905.1_KSU_PyrO_219_1","CM091585.1",1410270,1427258),"MAT1-2":("GCA_905067035.2_AG039_contigs_polished","CAJHII020000004.1",2806290,2820991)},
 "botrytis_cinerea":{"MAT1-1":("GCF_000143535.2_ASM14353v4","NC_037310.1",808684,821303),"MAT1-2":("GCA_015148905.1_ASM1514890v1","JACVGA010000339.1",15977,28871)},
 "verticillium_dahliae":{"MAT1-1":("GCA_018140825.1_ASM1814082v1","MU211481.1",3995601,4031304),"MAT1-2":("GCA_000400815.2_VDAG_JR2v.4.0","CP009079.1",3830361,3865909)},
 "aspergillus_flavus":{"MAT1-1":("GCF_009017415.1_ASM901741v1","NC_092410.1",2426714,2458292),"MAT1-2":("GCA_014117485.1_ASM1411748v1","CP059863.1",1617807,1649488)},
}
def contig(path, name):
    seq=[]; on=False
    with gzip.open(path,"rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on: break
                on = line[1:].split()[0]==name; continue
            if on: seq.append(line.strip())
    return "".join(seq)
for sp,spec in SPEC.items():
    od=OUT+sp+"/"; os.makedirs(od,exist_ok=True); pre=sp.split("_")[0]
    for idio,(g,c,s,e) in spec.items():
        sq=contig(D+g+".fa.gz",c); assert sq,(g,c)
        sub=sq[s-1:e]
        open(f"{od}{pre}_{idio}.fasta","w").write(f">{sp}_{idio} {g}:{c}:{s}-{e}\n"+"\n".join(sub[i:i+80] for i in range(0,len(sub),80))+"\n")
    p=Panel.from_fastas({i:f"{od}{pre}_{i}.fasta" for i in spec},k=31,min_run=100)
    print(sp,{i:p.lengths[i] for i in spec},"unique",{k:len(v) for k,v in p.unique.items()},"shared",len(p.shared),flush=True)
