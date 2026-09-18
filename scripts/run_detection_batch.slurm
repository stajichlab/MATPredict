#!/usr/bin/bash -l
#SBATCH -N 1 -c 4 --mem 16gb --time=0-02:30:00 --out logs/run_detection_batch.%A.log
# Runs one detection batch produced by MATPredict.detect.batch_runner.plan_batches.
#
# Usage:
#   sbatch scripts/run_detection_batch.slurm <genomes_json> <db_root> <out_dir>
#
# <genomes_json> is one batch's genome list (see scripts/run_detection_batch.py's
# docstring for the expected JSON shape). Each SLURM job is sized toward
# plan_batches's target_seconds_per_job (default 5400s / 1.5h); the 2h30m wall
# time above is headroom above that target, not the expected runtime.

: "${SCRATCH:?SCRATCH must be set by SLURM for this job}"

GENOMES_JSON="$1"
DB_ROOT="$2"
OUT_DIR="$3"

if [[ -z "$GENOMES_JSON" || -z "$DB_ROOT" || -z "$OUT_DIR" ]]; then
    echo "usage: sbatch run_detection_batch.slurm <genomes_json> <db_root> <out_dir>" >&2
    exit 1
fi

pixi run python scripts/run_detection_batch.py \
    --genomes-json "$GENOMES_JSON" \
    --db-root "$DB_ROOT" \
    --out-dir "$OUT_DIR"
