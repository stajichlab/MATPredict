import csv,sys
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates/synteny"
lab=sys.argv[1]
rows=[r for r in csv.DictReader(open(O+"/pjir_flank_orthologs.tsv"),delimiter="\t") if r["target"]==lab]
rows.sort(key=lambda r:int(r["pjir_start"]))
for r in rows: print(f'Pjir {int(r["pjir_start"])//1000:>4}kb {r["side"]:2s} {r["pjir_gene"][-10:]}  ->  {r["contig"]}:{int(r["start"])//1000}kb {r["strand"]} {float(r["pident"]):.0f}% e={r["evalue"]}')
