import gzip, sys, os
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-b317ece/src")
from MATPredict.reads.panel import Panel, read_fasta_seq
from MATPredict.reads.kmers import canonical_kmers
D="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/panels/pyricularia_oryzae/"
def contig(path, name):
    seq=[]; on=False
    with gzip.open(path,"rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on: break
                on = line[1:].split()[0]==name; continue
            if on: seq.append(line.strip())
    return "".join(seq)
extra={"MAT1-1":[("GCA_031208725.1_ASM3120872v1","JASWIH010000765.1",113323,150445)],
       "MAT1-2":[("GCA_046718435.1_UK_UFVPY232","JBCEVS010001119.1",38568,52859),("GCA_047323855.1_ASM4732385v1","JBDQYV010000213.1",6169,11294)]}
base={i:read_fasta_seq(O+f"pyricularia_{i}.fasta") for i in ("MAT1-1","MAT1-2")}
p=Panel.from_fastas({i:O+f"pyricularia_{i}.fasta" for i in base},k=31,min_run=100)
allp=p.shared|p.unique["MAT1-1"]|p.unique["MAT1-2"]
for idio,lst in extra.items():
    for g,c,s,e in lst:
        sq=contig(D+g+".fa.gz",c)[s-1:e]; ks=canonical_kmers(sq,31)
        print(f"{idio} {g} {c}:{s}-{e} len {len(sq)}: kmers in panel: shared {len(ks&p.shared)/max(1,len(p.shared)):.2f} (of panel shared), unique1 {len(ks&p.unique['MAT1-1'])/len(p.unique['MAT1-1']):.2f}, unique2 {len(ks&p.unique['MAT1-2'])/len(p.unique['MAT1-2']):.2f}; fraction of its kmers in panel {len(ks&allp)/len(ks):.2f}")
