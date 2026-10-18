#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J doth-scan -p short
# array task N scans scan_chunks/scan_NNN.tsv (ASM_ACCESSION column) with scan_genome.py -> scan/scan_NNN.tsv
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-c9c0ec6
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test/bin:$PATH
export PYTHONPATH=$REPO/src PYTHONNOUSERSITE=1
OUT=$REPO/results/2026-10-18_dothideomycetes_audit
cd $OUT; mkdir -p scan
N=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
[ -s scan/scan_$N.tsv ] && [ "$(wc -l < scan/scan_$N.tsv)" -gt "$(wc -l < scan_chunks/scan_$N.tsv)" ] && { echo done already; exit 0; }
python scan_genome.py --list scan_chunks/scan_$N.tsv --lib /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes \
   --db $REPO/db --pfam-subset pfam_subset.hmm --out scan/scan_$N.tsv
