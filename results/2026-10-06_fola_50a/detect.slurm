#!/bin/bash
#SBATCH -p exfab -A exfab -c 8 --mem 24gb -t 2:00:00 -J fola50a-detect --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_50a/detect_%j.log
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-94c1a3b/src MATPREDICT_DB_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-94c1a3b/db MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache
python -m MATPredict detect --genome /bigdata/stajichlab/jstajich/projects/fola_reads/50a_spades/scaffolds.fasta --taxid 5507 --out-dir /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_50a/detect
