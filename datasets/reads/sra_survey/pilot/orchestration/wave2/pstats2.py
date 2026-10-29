import sys
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-b317ece/src")
from MATPredict.reads.panel import Panel
B="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/panels_cf/"
for d,pre in (("pyricularia_oryzae","pyricularia"),("cryphonectria_parasitica","cryphonectria"),("zymoseptoria_tritici","zymoseptoria")):
    p=Panel.from_fastas({i:f"{B}{d}/{pre}_{i}.fasta" for i in ("MAT1-1","MAT1-2")},k=31,min_run=100)
    print(d,"lengths",p.lengths,"unique",{k:len(v) for k,v in p.unique.items()},"shared",len(p.shared),flush=True)
