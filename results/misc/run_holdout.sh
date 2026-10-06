#!/usr/bin/bash
set -u
W=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts
[[ -f "$W/src/MATPredict/__init__.py" ]] || { echo "SRC missing on $(hostname)" >&2; exit 1; }
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python \
  $W/scripts/run_holdout_benchmark.py \
  --worktree $W \
  --assemblies /bigdata/stajichlab/jstajich/projects/MATPredict/results/truth_assemblies.json \
  --work ${SCRATCH:?}/holdout \
  --out /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-23_holdout_benchmark.yaml
