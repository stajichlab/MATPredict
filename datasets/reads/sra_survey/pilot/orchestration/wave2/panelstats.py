import sys
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-b317ece/src")
from MATPredict.reads.panel import Panel
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/panel_ncrassa/"
for mr in (100,):
    p=Panel.from_fastas({"MAT1-1":O+"neurospora_MAT1-1.fasta","MAT1-2":O+"neurospora_MAT1-2.fasta"},k=31,min_run=mr)
    print("min_run",mr,"lengths",p.lengths,"unique",{k:len(v) for k,v in p.unique.items()},"shared",len(p.shared))
