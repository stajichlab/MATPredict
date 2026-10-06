import gzip,sys,os
from multiprocessing import Pool
D="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
def st(a):
    p=D+a+".fa.gz"
    if not os.path.exists(p): return f"{a}\tNA\tNA\tNA\tNA"
    L=[];cur=0
    with gzip.open(p,'rt') as f:
        for line in f:
            if line[0]=='>':
                if cur: L.append(cur)
                cur=0
            else: cur+=len(line.strip())
    if cur: L.append(cur)
    L.sort(reverse=True); tot=sum(L); c=0
    for x in L:
        c+=x
        if c>=tot/2: n50=x;break
    return f"{a}\t{len(L)}\t{tot}\t{n50}\t{L[0]}"
ids=[l.strip() for l in open("asmids.txt")]
with Pool(4) as p:
    with open("asmstats.tsv","w") as o:
        o.write("asm\tn_contigs\ttotal\tN50\tmax\n")
        for r in p.imap(st,ids,chunksize=8): o.write(r+"\n")
