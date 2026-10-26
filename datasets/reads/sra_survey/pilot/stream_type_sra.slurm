#!/bin/bash -l
#SBATCH -p exfab -A exfab -c 4 --mem 8gb -t 2:00:00 -J sra-type
# Slurm wrapper: one array task = one row of RUNS.tsv, run through genotype_mat_one_accession.sh (which knows nothing about Slurm; a Nextflow process can call it the same way).
# Usage: sbatch --array=1-N%8 stream_type_sra.slurm RUNS.tsv PANEL_DIR PANEL_PREFIX OUT_DIR
#   RUNS.tsv: tab-separated with a header (see CONTRACT.md): run_accession, read_count; optional group, library_layout, panel_dir, panel_prefix.
#   panel_dir / panel_prefix in a row override the two arguments, so one array can mix species.
# Env: TARGET_SPOTS, WINDOWS, SRA_URL, KEEP_SRA, IDIOMORPHS are passed on to genotype_mat_one_accession.sh.
set -uo pipefail
RUNS=$1; PDIR=$2; PFX=$3; OUT=$4
here=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" 2>/dev/null && pwd)
[[ -x $here/genotype_mat_one_accession.sh ]] || here=${WORKER_DIR:-/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit}
: "${SCRATCH:?SCRATCH must be set by Slurm}"; export SCRATCH
getf() { awk -F'\t' -v name="$1" -v row=$((SLURM_ARRAY_TASK_ID+1)) 'NR==1{for(i=1;i<=NF;i++)if($i==name)c=i} NR==row{print (c?$c:"")}' "$RUNS"; }
run=$(getf run_accession); [[ -n $run ]] || { echo "no run for task $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
rp=$(getf panel_dir); rx=$(getf panel_prefix); [[ -n $rp ]] && PDIR=$rp; [[ -n $rx ]] && PFX=$rx
exec "$here/genotype_mat_one_accession.sh" "$run" "$(getf read_count)" "$PDIR" "$PFX" "$OUT" "$(getf group)" "$(getf library_layout)"
