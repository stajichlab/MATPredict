#!/usr/bin/env bash
# download.sh - fetch one ENA run (both mates) per Slurm array task, verify md5.
#
# Usage (Slurm array, one run per task; N = number of data rows in selected.tsv):
#   N=$(( $(wc -l < selected.tsv) - 1 ))
#   sbatch --array=1-${N}%8 -c1 --mem=2G -t 06:00:00 -J crypto_reads \
#          --wrap "bash $PWD/download.sh /path/to/outdir $PWD/selected.tsv"
# or put in a batch file:   #SBATCH --array=1-312%8
#   bash download.sh OUTDIR [TSV]
#
# Args:  $1 = target directory (required; files go to $1/<run_accession>/)
#        $2 = TSV with header incl. run_accession, fastq_ftp, fastq_md5 (default: selected.tsv beside this script)
# Task selection: SLURM_ARRAY_TASK_ID (1-based data row). Outside Slurm, set TASK_ID=<n>,
#   or TASK_ID=all to loop over every row serially.
# Idempotent: a run is skipped if $OUT/<run>/.done exists; partial files are resumed (curl -C -);
#   a file whose md5 already matches is not re-downloaded. Keep %N concurrency small to be polite to ENA.
# Env: ENA_PROTO=https (default; ftp also works), RETRIES=5.
set -euo pipefail

OUT=${1:?usage: download.sh OUTDIR [TSV]}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TSV=${2:-$HERE/selected.tsv}
PROTO=${ENA_PROTO:-https}
RETRIES=${RETRIES:-5}
TASK=${SLURM_ARRAY_TASK_ID:-${TASK_ID:-}}
[[ -n $TASK ]] || { echo "No SLURM_ARRAY_TASK_ID or TASK_ID set" >&2; exit 2; }
[[ -r $TSV ]] || { echo "Cannot read $TSV" >&2; exit 2; }
mkdir -p "$OUT"

if command -v md5sum >/dev/null 2>&1; then md5of() { md5sum "$1" | cut -d' ' -f1; }
else md5of() { md5 -q "$1"; }; fi

fetch_row() {  # $1 = data row number (1-based)
  local line run ftps md5s
  line=$(awk -F'\t' -v n="$1" '
    NR==1 {for(i=1;i<=NF;i++){if($i=="run_accession")r=i; if($i=="fastq_ftp")f=i; if($i=="fastq_md5")m=i}
           if(!r||!f||!m){print "ERR missing columns" > "/dev/stderr"; exit 3}; next}
    NR==n+1 {print $r "\t" $f "\t" $m; exit}' "$TSV") || return 3
  [[ -n $line ]] || { echo "No row $1 in $TSV" >&2; return 2; }
  IFS=$'\t' read -r run ftps md5s <<<"$line"
  local d="$OUT/$run"
  if [[ -e $d/.done ]]; then echo "[$run] already done"; return 0; fi
  mkdir -p "$d"
  IFS=';' read -r -a U <<<"$ftps"; IFS=';' read -r -a M <<<"$md5s"
  [[ ${#U[@]} -eq ${#M[@]} && ${#U[@]} -ge 1 ]] || { echo "[$run] url/md5 count mismatch" >&2; return 4; }
  local i f want got
  for i in "${!U[@]}"; do
    f="$d/$(basename "${U[$i]}")"; want=${M[$i]}
    if [[ -s $f ]] && [[ $(md5of "$f") == "$want" ]]; then echo "[$run] $(basename "$f") ok (cached)"; continue; fi
    echo "[$run] downloading $(basename "$f")"
    curl -fL --retry "$RETRIES" --retry-delay 10 --retry-connrefused -C - -o "$f" "$PROTO://${U[$i]}" \
      || { echo "[$run] curl failed (partial kept for resume)" >&2; return 5; }
    got=$(md5of "$f")
    if [[ $got != "$want" ]]; then
      echo "[$run] MD5 MISMATCH $(basename "$f"): got $got want $want" >&2
      rm -f "$f"; return 6
    fi
    echo "[$run] $(basename "$f") md5 ok"
  done
  date -u +%FT%TZ > "$d/.done"
}

if [[ $TASK == all ]]; then
  n=$(( $(wc -l < "$TSV") - 1 )); rc=0
  for ((k=1;k<=n;k++)); do fetch_row "$k" || rc=1; done
  exit $rc
else
  fetch_row "$TASK"
fi
