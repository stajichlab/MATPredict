import sys
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-b317ece/src")
from MATPredict.reads.panel import Panel
B="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/panels_cf/sclerotinia_sclerotiorum/"
p=Panel.from_fastas({i:f"{B}sclerotinia_{i}.fasta" for i in ("MAT1-1","MAT1-2")},k=31,min_run=100)
print("lengths",p.lengths,"unique",{k:len(v) for k,v in p.unique.items()},"shared",len(p.shared))
