import csv,os
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave1"
def rd(p): return list(csv.DictReader(open(p,newline="").read().replace("\r","").splitlines(),delimiter="\t"))
done={f.split(".")[0] for f in os.listdir(W+"/out") if f.endswith(".run.json")}
rows=rd(W+"/plan/clonostachys_rosea/runs_full.tsv")[:50]; todo=[r for r in rows if r["run_accession"] not in done]
w=csv.DictWriter(open("full3_runs.tsv","w"),list(rows[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(todo); print(len(todo))
