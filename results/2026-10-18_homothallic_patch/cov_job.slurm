#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 01:00:00 -p short -J hl-cov
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/hl-patch
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test/bin:$PATH PYTHONPATH=$REPO/src PYTHONNOUSERSITE=1
cd $REPO/results/2026-10-18_homothallic_patch
python cov_check.py $REPO/db runs_out/runs /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes cov_check.tsv
