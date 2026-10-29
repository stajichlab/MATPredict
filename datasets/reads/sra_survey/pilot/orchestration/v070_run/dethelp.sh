#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict
export PYTHONPATH=$P/.claude/worktrees/run-b317ece/src MATPREDICT_DB_ROOT=$P/.claude/worktrees/run-b317ece/db
$P/.pixi/envs/default/bin/python -m MATPredict detect --help 2>&1 | grep -E "^\s+--(taxid|phylum|family|order|genome|out|preset|exhaustive|lineage|db|records|exclude|min)" | cut -c1-170
du -sh $P/.claude/worktrees/run-b317ece/db | cut -f1
grep -c . $P/.claude/worktrees/run-b317ece/db/Ascomycota/Taphrinales/5011_unknown-1_PM_combined/proteins.faa
head -4 $P/.claude/worktrees/run-b317ece/db/Ascomycota/Taphrinales/5011_unknown-1_PM_combined/proteins.faa | cut -c1-120
