#!/bin/bash
#SBATCH -p exfab -A exfab -c 3 --mem 12gb -t 12:00:00 -J afum-detect --array=0-9 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/detect/slurm_%A_%a.log
# matpredict detect (frozen run-2d1860a = main at PR #41) on all 304 scaffolded A. fumigatus genomes; species taxid 746128.
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-2d1860a/src MATPREDICT_DB_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-2d1860a/db MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache
run(){ g=$1; s=$(basename $g .sorted.fasta); python -m MATPredict detect --genome $g --taxid 746128 --out-dir /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/detect/runs/$s > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/detect/runs/$s.log 2>&1 || echo "FAIL $s"; }
export -f run
awk -v t=$SLURM_ARRAY_TASK_ID 'NR%10==t' /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/genomes.txt | xargs -P3 -I{} bash -c 'run {}'
echo done
