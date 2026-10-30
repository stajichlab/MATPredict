import csv
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave1"
slugs=["histoplasma_capsulatum","neurospora_tetraspora","aspergillus_latus","fusarium_avenaceum","lasallia_pustulata","exserohilum_rostratum","trichoderma_koningiopsis"]
rows=[]
for s in slugs:
    rows+=list(csv.DictReader(open(f"{R}/plan/{s}/runs_full.tsv"),delimiter="\t"))[:50]
w=csv.DictWriter(open(f"{R}/runs_full2.tsv","w"),rows[0].keys(),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
print(len(rows),"runs for",len(slugs),"species")
