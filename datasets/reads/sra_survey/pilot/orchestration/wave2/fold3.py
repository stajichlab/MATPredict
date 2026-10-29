import csv,sys,shutil
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave1"
rd=lambda p:list(csv.DictReader(open(p,newline="").read().replace("\r","").splitlines(),delimiter="\t"))
fin=rd(W+"/final/summary/results_all.tsv"); have={r["run"] for r in fin}
new=[r for r in rd(W+"/full3/summary/results_all.tsv") if r["run"] not in have]
cols=list(fin[0].keys())
for r in new:
    for c in cols: r.setdefault(c,"")
shutil.copy(W+"/final/summary/results_all.tsv",W+"/final/summary/results_all.before_full3.tsv")
with open(W+"/final/summary/results_all.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,cols,delimiter="\t",lineterminator="\n",extrasaction="ignore"); w.writeheader(); w.writerows(fin+new)
print("added",len(new),"runs; total",len(fin)+len(new))
