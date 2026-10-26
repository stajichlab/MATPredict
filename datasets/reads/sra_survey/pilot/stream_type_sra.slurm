#!/bin/bash -l
#SBATCH -p exfab -A exfab -c 4 --mem 8gb -t 2:00:00 -J sra-type
# One SRA run per array task: fetch the run's .sra from NCBI's copy on AWS open data (no login, no cost to us), cut evenly spaced
# windows of spots with fastq-dump (random access), type them against a species panel, keep only the small result files, delete the reads.
# Why AWS and not ENA: from the HPCC (2026-10-26) one ENA stream ran at 0.1-0.9 MB/s, the AWS copy at about 20 MB/s per stream
# (76 MB/s with four). The windows also fix position-sorted old runs: reads are taken along the whole run, not from its start.
#
# Usage: sbatch --array=1-20%4 stream_type_sra.slurm RUNS.tsv PANEL_DIR PANEL_PREFIX OUT_DIR
#   RUNS.tsv: tab-separated with a header; needs run_accession and read_count (spots), optional group, library_layout, study_accession.
#   Panel files: PANEL_DIR/PANEL_PREFIX_<IDIOMORPH>.fasta; IDIOMORPHS="MAT1-1 MAT1-2" (default) names the idiomorphs.
# Env: TARGET_SPOTS (default 2000000 spots, 2 reads each for paired runs), WINDOWS (8), SRA_URL, KEEP_SRA=1 to keep the .sra.
# Kept per run in OUT_DIR: <run>.<prefix>.tsv (reads-type row), <run>.run.json (accession, layout, spots, windows, bytes, seconds, state).
# Not kept: the .sra and every FASTQ (removed by the exit trap).
set -uo pipefail
RUNS=$1; PDIR=$2; PFX=$3; OUT=$4
TARGET=${TARGET_SPOTS:-2000000}; WIN=${WINDOWS:-8}
SRA_URL=${SRA_URL:-https://sra-pub-run-odp.s3.amazonaws.com/sra}
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/sratoolkit/3.2.0/bin:$PATH
P=/bigdata/stajichlab/jstajich/projects/MATPredict
PY=$P/.pixi/envs/default/bin/python
export PYTHONPATH=$P/.claude/worktrees/run-b317ece/src
export MATPREDICT_DB_ROOT=$P/.claude/worktrees/run-b317ece/db
: "${SCRATCH:?SCRATCH must be set by Slurm}"
mkdir -p "$OUT"
getf() { awk -F'\t' -v name="$1" -v row=$((SLURM_ARRAY_TASK_ID+1)) 'NR==1{for(i=1;i<=NF;i++)if($i==name)c=i} NR==row{print (c?$c:"")}' "$RUNS"; }
run=$(getf run_accession); spots=$(getf read_count); layout=$(getf library_layout); group=$(getf group)
[[ -n $run ]] || { echo "no run for task $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
[[ -s $OUT/$run.$PFX.tsv ]] && { echo "skip $run"; exit 0; }
work=$(mktemp -d "$SCRATCH/sra_type.$run.XXXXXX"); trap '[[ ${KEEP_SRA:-0} == 1 ]] || rm -rf "$work"' EXIT
state=started; t0=$(date +%s)
json() { printf '{"run":"%s","group":"%s","layout":"%s","spots_total":%s,"target_spots":%s,"windows":%s,"sra_bytes":%s,"secs_download":%s,"secs_convert":%s,"secs_type":%s,"state":"%s","source":"aws-open-data","host":"%s"}\n' \
  "$run" "$group" "$layout" "${spots:-0}" "$TARGET" "$WIN" "${bytes:-0}" "${d1:-0}" "${c1:-0}" "${y1:-0}" "$state" "$(hostname)" > "$OUT/$run.run.json"; }
trap 'json; [[ ${KEEP_SRA:-0} == 1 ]] || rm -rf "$work"' EXIT
# 1. download (resumable, 3 tries)
ok=0; for a in 1 2 3; do
  curl -sS --fail -L -C - -o "$work/$run.sra" "$SRA_URL/$run/$run" && { ok=1; break; }; sleep $((a*10)); done
d1=$(( $(date +%s) - t0 )); bytes=$(stat -c %s "$work/$run.sra" 2>/dev/null || echo 0)
[[ $ok == 1 && $bytes -gt 0 ]] || { state=failed_download; json; exit 1; }
# 2. windows of spots (1-based, inclusive); a run shorter than the target is read whole
t1=$(date +%s); fq="$work/sub.fastq"; : > "$fq"
if [[ -z $spots || $spots -le $TARGET ]]; then
  fastq-dump -Z --split-spot --ncbi_error_report never "$work/$run.sra" 2>/dev/null >> "$fq"
else
  w=$(( TARGET / WIN ))
  for i in $(seq 0 $((WIN-1))); do
    s=$(( (2*i+1) * spots / (2*WIN) - w/2 + 1 )); [[ $s -lt 1 ]] && s=1; e=$(( s + w - 1 ))
    fastq-dump -N $s -X $e -Z --split-spot --ncbi_error_report never "$work/$run.sra" 2>/dev/null >> "$fq"
  done
fi
c1=$(( $(date +%s) - t1 )); nlines=$(wc -l < "$fq")
[[ $nlines -gt 0 ]] || { state=failed_convert; json; exit 1; }
rm -f "$work/$run.sra"   # the archive is not needed any more; only the subsample stays until the trap
# 3. type
t2=$(date +%s)
printf 'sample\treads\n%s\t%s\n' "$run" "$fq" > "$work/samples.tsv"
idio=(); for n in ${IDIOMORPHS:-MAT1-1 MAT1-2}; do idio+=(--idiomorph "$n=$PDIR/${PFX}_$n.fasta"); done
$PY -m MATPredict reads-type "${idio[@]}" --samples "$work/samples.tsv" --out "$work/$run.$PFX.tsv" || { state=failed_type; y1=$(( $(date +%s) - t2 )); json; exit 1; }
y1=$(( $(date +%s) - t2 )); state=done
mv "$work/$run.$PFX.tsv" "$OUT/$run.$PFX.tsv"
echo "done $run: $((nlines/4)) reads typed, download ${d1}s ($((bytes/1000000)) MB), convert ${c1}s, type ${y1}s"
