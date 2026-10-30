import csv
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates/synteny"
R=list(csv.DictReader(open(O+"/pjir_flank_orthologs.tsv"),delimiter="\t"))
for lab in ("Pcar_B80","Tdef_PYCC","Scoloradoensis","Scomp","Neolecta","Plactucae","Spombe"):
    print("==",lab)
    for r in sorted([x for x in R if x["target"]==lab and 200000<=int(x["pjir_start"])<=236000],key=lambda x:int(x["pjir_start"])):
        print(f'   Pjir {int(r["pjir_start"])//1000:>4}kb {r["side"]:2s} -> {r["contig"]}:{int(r["start"])//1000}kb{r["strand"]} {float(r["pident"]):.0f}% e={r["evalue"]}')
