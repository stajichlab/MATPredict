#!/bin/bash
#SBATCH -p exfab -A exfab -c 16 --mem 48gb -t 4:00:00 -J fola-detect2 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_detect_all/slurm_%j.log
# MATPredict detect, frozen run-94c1a3b (F. oxysporum MAT1-1/MAT1-2 records), taxid 5507, all Fola assemblies + VSP-0947.
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-94c1a3b/src MATPREDICT_DB_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-94c1a3b/db MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache
run(){ f=$1; s=${f%.*}; python -m MATPredict detect --genome /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_detect_all/genomes/$f --taxid 5507 --out-dir /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_detect_all/runs/$s > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_detect_all/runs_$s.log 2>&1 || echo "FAIL $s"; }
export -f run
ls /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_detect_all/genomes | xargs -P 16 -I{} bash -c 'run {}'
echo done
