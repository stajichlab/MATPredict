import csv,sys
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave1/plan"
out=csv.writer(open(sys.argv[1],"w"),delimiter="\t",lineterminator="\n"); n=0
hdr=None
for sp in ["fusarium_oxysporum","aspergillus_fumigatus","trichoderma_reesei","penicillium_roqueforti","coccidioides_posadasii","verticillium_albo_atrum"]:
    for i,r in enumerate(csv.DictReader(open(f"{W}/{sp}/runs_full.tsv"),delimiter="\t")):
        if hdr is None: hdr=list(r); out.writerow(hdr)
        out.writerow([r[h] for h in hdr]); n+=1
print(n)
