#!/bin/bash
#SBATCH -p exfab -A exfab -c 2 --mem 8gb -t 3:00:00 -J mp-tests-ss2 --array=0-7 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_supported_span_tests2/out/slurm_%A_%a.log
set +eu; source /etc/profile; set -eu
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/core-span/src PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/.pixi/envs/test/bin:$PATH
cd /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/core-span
mapfile -t F < <(awk -v t=$SLURM_ARRAY_TASK_ID 'NR%8==t' /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_supported_span_tests2/test_files.txt)
python -m pytest "${F[@]}" -q -p no:cacheprovider
echo "shard $SLURM_ARRAY_TASK_ID exit $?"
