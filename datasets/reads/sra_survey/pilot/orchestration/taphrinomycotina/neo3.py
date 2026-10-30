import csv
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates/loose"
for x in csv.DictReader(open(O+"/candidates.tsv"),delimiter="\t"):
    if x["species"].startswith(("Neolecta","Schizosaccharomyces japonicus","Saitoella")): print({k:x[k] for k in ("species","contig","start","end","best_ref","ref_class","identity","query_coverage","protein_aa","report_locus")})
