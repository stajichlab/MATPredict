#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict; DB=$P/.claude/worktrees/run-b317ece/db/Ascomycota
export PATH=$P/.pixi/envs/default/bin:$PATH
T=$(mktemp -d /tmp/dc.XXXX)
cat $DB/Dipodascales/4952_clib122_MATyl_B/proteins.faa $DB/Dipodascales/4952_w29_MATyl_A/proteins.faa > $T/yl.faa 2>/dev/null
ls $DB/Dipodascales | head -8 | tr '\n' ' '; echo
for d in $DB/*/1173061_*; do cat $d/proteins.faa; done > $T/gc.faa
grep ">" $T/yl.faa | sed 's/.*name=//' | tr '\n' ' '; echo; grep ">" $T/gc.faa | sed 's/.*name=//' | tr '\n' ' '; echo
makeblastdb -in $T/yl.faa -dbtype prot -out $T/yl -logfile /dev/null
blastp -query $T/gc.faa -db $T/yl -evalue 1e-3 -outfmt "6 qseqid sseqid pident length evalue bitscore" -max_target_seqs 6 | sed 's/|gene_index=[0-9]*|name=/:/g;s/|role=[a-zA-Z_]*//g' | sort -k1,1 -k6,6nr | head -14
