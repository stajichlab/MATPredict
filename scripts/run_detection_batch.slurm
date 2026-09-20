#!/usr/bin/bash -l
#SBATCH -N 1 -c 4 --mem 16gb --time=0-02:30:00 --out logs/run_detection_batch.%A.log
# Runs one detection batch produced by MATPredict.detect.batch_runner.plan_batches.
#
# Usage:
#   sbatch scripts/run_detection_batch.slurm <genomes_json> <db_root> <out_dir> [driver flags...]
#
# Examples:
#   sbatch scripts/run_detection_batch.slurm batch_003.json db out/batch_003
#   sbatch scripts/run_detection_batch.slurm batch_003.json db out/batch_003 \
#       --phylum Mucoromycota
#   sbatch scripts/run_detection_batch.slurm batch_003.json db out/batch_003 \
#       --phylum Mucoromycota --min-hits 1 --no-require-core-role
#
# <genomes_json> is one batch's genome list (see scripts/run_detection_batch.py's
# docstring for the expected JSON shape). Each SLURM job is sized toward
# plan_batches's target_seconds_per_job (default 5400s / 1.5h); the 2h30m wall
# time above is headroom above that target, not the expected runtime.
#
# Everything after the three positionals is forwarded VERBATIM to
# scripts/run_detection_batch.py, so every driver flag works here with the same
# spelling it has in `matpredict detect`: --phylum, --min-hits, --min-identity,
# --require-core-role / --no-require-core-role. Forwarding the tail rather than
# re-declaring each flag is deliberate -- this wrapper is the documented rollout
# entry point, and a hand-maintained flag list here would silently drop any flag
# the driver gains, which is exactly how an 813-genome rollout ends up searching
# all 181 curated proteins per genome after the narrowing work was done.

: "${SCRATCH:?SCRATCH must be set by SLURM for this job}"

GENOMES_JSON="$1"
DB_ROOT="$2"
OUT_DIR="$3"

if [[ -z "$GENOMES_JSON" || -z "$DB_ROOT" || -z "$OUT_DIR" ]]; then
    echo "usage: sbatch run_detection_batch.slurm <genomes_json> <db_root> <out_dir> [driver flags...]" >&2
    exit 1
fi

shift 3
DRIVER_ARGS=("$@")

# The scope this job actually ran with is printed into the job log. Whether the
# query set was narrowed changes what the run's negatives mean, and an operator
# reading a finished job's log must be able to tell WITHOUT re-deriving it from
# the sbatch command line, which is not recorded anywhere in the log.
SCOPED_PHYLUM=""
for ((i = 0; i < ${#DRIVER_ARGS[@]}; i++)); do
    case "${DRIVER_ARGS[i]}" in
        --phylum) SCOPED_PHYLUM="${DRIVER_ARGS[i + 1]}" ;;
        --phylum=*) SCOPED_PHYLUM="${DRIVER_ARGS[i]#--phylum=}" ;;
    esac
done

if [[ -n "$SCOPED_PHYLUM" ]]; then
    echo "scope: phylum=${SCOPED_PHYLUM} (routing and the tblastn query set are restricted to this phylum)"
else
    echo "scope: phylum=<none> (per-genome self-routing; the query set is every curated protein)"
fi
echo "driver flags: ${DRIVER_ARGS[*]:-<none>}"

pixi run python scripts/run_detection_batch.py \
    --genomes-json "$GENOMES_JSON" \
    --db-root "$DB_ROOT" \
    --out-dir "$OUT_DIR" \
    "${DRIVER_ARGS[@]}"
