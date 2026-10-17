#!/bin/bash
#SBATCH -p short -A stajichlab -c 4 --mem 8gb -t 4:00:00 -J fetch-type
# One SRA run per array task: download to $SCRATCH, type against a species panel, keep only the result row, delete the reads.
# Nothing but the small result TSV lands under /bigdata.
#
# Usage:
#   N=$(( $(wc -l < selected.tsv) - 1 ))
#   sbatch --array=1-${N}%8 fetch_type.slurm <selected.tsv> <downloader.sh> <panel_dir> <panel_prefix> <out_dir>
# e.g. sbatch --array=1-312%8 fetch_type.slurm cryptococcus/selected.tsv cryptococcus/download.sh panels/cryptococcus neoformans out/crypto
# <panel_dir>/<panel_prefix>_MAT1-1.fasta and _MAT1-2.fasta must exist (add more --idiomorph panels by editing the call).
# Env: source datasets/reads/coccidioides/pilot/env.sh first (or set CLONE/TOOLS the same way).
set -uo pipefail
: "${SCRATCH:?SCRATCH must be set by Slurm}"
tsv=$1; dl=$2; pdir=$3; pfx=$4; out=$5
run=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1}' "$tsv")
mkdir -p "$out"
[[ -s $out/$run.tsv ]] && { echo "skip $run"; exit 0; }
work=$(mktemp -d "$SCRATCH/fetch_type.$run.XXXXXX"); trap 'rm -rf "$work"' EXIT
bash "$dl" "$tsv" "$work/reads" "$run" || { echo "FAIL download $run" >&2; exit 1; }
printf 'sample\treads\n%s\t%s,%s\n' "$run" "$work/reads/$run/${run}_1.fastq.gz" "$work/reads/$run/${run}_2.fastq.gz" > "$work/samples.tsv"
python -m MATPredict reads-type --idiomorph MAT1-1="$pdir/${pfx}_MAT1-1.fasta" --idiomorph MAT1-2="$pdir/${pfx}_MAT1-2.fasta" \
  --samples "$work/samples.tsv" --max-reads 8000000 --out "$work/$run.tsv" || { echo "FAIL type $run" >&2; exit 1; }
mv "$work/$run.tsv" "$out/$run.tsv"   # result only; the trap removes the reads
