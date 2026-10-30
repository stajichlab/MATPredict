import csv
P="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2/botrytis/runs_calls.tsv"
r=list(csv.DictReader(open(P),delimiter="\t"))
print([ (x["run"],x["strain"],x["study"],x["call"],x["truth"],x["shared_depth"]) for x in r if x["truth"] in ("MAT1-1","MAT1-2") and x["call"] in ("MAT1-1","MAT1-2") and x["call"]!=x["truth"]])
print("with truth:",sum(1 for x in r if x["truth"]),"of",len(r))
from math import comb
n,k=33,18; p=sum(comb(n,i) for i in range(0,n+1) if abs(i-n/2)>=abs(k-n/2))/2**n; print("two-sided binomial p vs 1:1 for 18 of 33: %.2f"%p)
