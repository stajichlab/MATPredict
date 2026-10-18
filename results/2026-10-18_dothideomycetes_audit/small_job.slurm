#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 01:30:00 -p short
# generic: small_job.slurm SCRIPT [args...]  -- runs a python script with the audit worktree's src and the `test` env
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-c9c0ec6
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test/bin:$PATH
export PYTHONPATH=$REPO/src PYTHONNOUSERSITE=1
export MATPREDICT_DB_ROOT=$REPO/db MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache
cd $REPO/results/2026-10-18_dothideomycetes_audit
python "$@"
