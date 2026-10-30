import csv,os,sys
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2"
rd=lambda p:list(csv.DictReader(open(p,newline="").read().replace("\r","").splitlines(),delimiter="\t"))
test=rd(W+"/test/summary/results_all.tsv"); full=rd(W+"/full/summary/results_all.tsv")
over={}
bp=W+"/botrytis_full/summary/results_all.tsv"
if os.path.exists(bp):
    for r in rd(bp): over.setdefault(r["group"],[]).append(r)
fgroups={r["group"] for r in full}
rows=[r for r in test if r["group"] not in fgroups and r["group"] not in over]+[r for r in full if r["group"] not in over]
for g,v in over.items(): rows+=v
cols=list(full[0].keys())
os.makedirs(W+"/final/summary",exist_ok=True)
with open(W+"/final/summary/results_all.tsv","w",newline="") as fh:
    w=csv.DictWriter(fh,cols,delimiter="\t",lineterminator="\n",extrasaction="ignore",restval=""); w.writeheader(); w.writerows(rows)
print("final rows",len(rows),"species",len({r["group"] for r in rows}),"override species:",{g:len(v) for g,v in over.items()})
