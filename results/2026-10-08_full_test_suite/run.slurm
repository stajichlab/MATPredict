#!/bin/bash
#SBATCH -p exfab -A exfab -c 2 --mem 8gb -t 3:00:00 -J mp-tests --array=0-7 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_full_test_suite/out/slurm_%A_%a.log
# Full pytest suite of main e9060da in 8 shards (test files assigned round-robin), frozen worktree + its pinned test env.
set +eu; source /etc/profile; set -eu
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/src
cd /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da
mapfile -t F < <(awk -v t=$SLURM_ARRAY_TASK_ID 'NR%8==t' /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_full_test_suite/test_files.txt)
echo "shard $SLURM_ARRAY_TASK_ID: ${#F[@]} files on $(hostname)"
start=$(date +%s)
/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/.pixi/envs/test/bin/python -m pytest "${F[@]}" -q -p no:cacheprovider --durations=15 -o junit_family=xunit2 --junitxml=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_full_test_suite/out/junit_$SLURM_ARRAY_TASK_ID.xml
rc=$?
echo "shard $SLURM_ARRAY_TASK_ID exit $rc wall $(($(date +%s)-start)) s"
