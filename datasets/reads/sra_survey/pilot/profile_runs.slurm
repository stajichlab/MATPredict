#!/bin/bash -l
#SBATCH -p exfab -A exfab -c 4 --mem 12gb -t 2:00:00 -J sra-profile
# Depth profile of one run along a panel (follow-up of a surprising call). Usage: sbatch --array=1-N profile_runs.slurm PROFILE_RUNS.tsv OUT_DIR
# PROFILE_RUNS.tsv columns: run_accession, read_count, panel_dir, panel_prefix, label. TARGET_SPOTS defaults to 6,000,000 here (more reads than a typing run).
set -uo pipefail
RUNS=$1; OUT=$2; here=${WORKER_DIR:-/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit}
getf() { awk -F'\t' -v name="$1" -v row=$((SLURM_ARRAY_TASK_ID+1)) 'NR==1{for(i=1;i<=NF;i++)if($i==name)c=i} NR==row{print (c?$c:"")}' "$RUNS"; }
run=$(getf run_accession); mkdir -p "$OUT"
fq="$SCRATCH/profile.$run.fastq"; trap 'rm -f "$fq"' EXIT
TARGET_SPOTS=${TARGET_SPOTS:-6000000} SKIP_TYPE=1 FASTQ_OUT=$fq SCRATCH=$SCRATCH "$here/genotype_mat_one_accession.sh" "$run" "$(getf read_count)" "$(getf panel_dir)" "$(getf panel_prefix)" "$OUT" profile || exit 1
P=/bigdata/stajichlab/jstajich/projects/MATPredict
PYTHONPATH=$P/.claude/worktrees/run-b317ece/src $P/.pixi/envs/default/bin/python "$here/depth_profile.py" --reads "$fq" --panel-dir "$(getf panel_dir)" --prefix "$(getf panel_prefix)" --out "$OUT/$run.$(getf label).profile.tsv"
