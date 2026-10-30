import csv,os,sys
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2"
rd=lambda p:list(csv.DictReader(open(p,newline="").read().replace("\r","").splitlines(),delimiter="\t"))
test=rd(W+"/test/summary/results_all.tsv"); full=rd(W+"/full/summary/results_all.tsv")
fullruns={r["run"] for r in full}; fgroups={r["group"] for r in full}
rows=[r for r in test if r["group"] not in fgroups]+full
cols=list(full[0].keys())
os.makedirs(W+"/final/summary",exist_ok=True)
with open(W+"/final/summary/results_all.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,cols,delimiter="\t",lineterminator="\n",extrasaction="ignore",restval=""); w.writeheader(); w.writerows(rows)
print("final rows",len(rows),"species",len({r["group"] for r in rows}),"from full:",len(full),"from test only:",len(rows)-len(full))
