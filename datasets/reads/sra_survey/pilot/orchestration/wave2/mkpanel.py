import gzip, sys, os
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-b317ece/src")
D="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
OUT="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/panel_ncrassa"
os.makedirs(OUT, exist_ok=True)
def contig(path, name):
    seq=[]; on=False
    with gzip.open(path,"rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on: break
                on = line[1:].split()[0]==name; continue
            if on: seq.append(line.strip())
    return "".join(seq)
regions={"MAT1-1":("GCF_000182925.2_NC12.fa.gz","NC_026501.1",1851649,1872342),
         "MAT1-2":("GCA_000786625.1_Neucr_trp3_1.fa.gz","KN389756.1",58545,77356)}
for name,(f,c,s,e) in regions.items():
    sq=contig(D+f,c); assert sq, c
    sub=sq[s-1:e]
    open(f"{OUT}/neurospora_{name}.fasta","w").write(f">neurospora_crassa_{name} {f.split('.fa')[0]}:{c}:{s}-{e}\n"+"\n".join(sub[i:i+80] for i in range(0,len(sub),80))+"\n")
    print(name, len(sub), "bp; contig", len(sq))
from MATPredict.reads import panel as P
print([n for n in dir(P) if not n.startswith("_")][:20])
