#!/bin/bash
#SBATCH -p exfab -A exfab -c 2 --mem 8gb -t 1:00:00 -J mp-tests-clf --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_full_test_suite/out/slurm_rerun_classifier_%j.log
# Re-run of the one file that failed in shard 2 only because mafft was not on PATH (the pixi env bin was not on PATH).
set +eu; source /etc/profile; set -eu
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/src PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/.pixi/envs/test/bin:$PATH
cd /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da
which mafft
python -m pytest tests/detect/test_classifier_build_deterministic.py -q -p no:cacheprovider --durations=5
