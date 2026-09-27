#!/usr/bin/bash
#SBATCH --partition=short --time=2:00:00 --cpus-per-task=8 --mem=16gb
#SBATCH --job-name=zygo23 --output=logs/zygo23_%j.log
# The Zygo 23 regression on both inputs (scaffold+proteome, contigs genome-only).
# Required: SRC=<frozen worktree>/src  OUT=<results dir>.  The db comes from the
# same tree (SRC/../db). Paths are passed in; nothing is resolved from this
# script's own location (that breaks under SLURM).
#   sbatch --export=ALL,SRC=$PWD/.claude/worktrees/run-<sha>/src,OUT=$PWD/results/<date>_zygo23 \
#       .claude/worktrees/run-<sha>/scripts/run_zygo_regression.slurm
set -euo pipefail
: "${SRC:?set SRC to a frozen worktree's src}" "${OUT:?set OUT}" "${SCRATCH:?no SCRATCH}"
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default
export PATH="$E/bin:$PATH"
[[ -f "$SRC/MATPredict/__init__.py" ]] || { echo "SRC missing: $SRC" >&2; exit 1; }
"$E/bin/python" "$(dirname "$SRC")/scripts/zygo_regression.py" run \
    --src "$SRC" --work "$SCRATCH/zygo23" --out "$OUT" --jobs "${SLURM_CPUS_PER_TASK:-8}"
