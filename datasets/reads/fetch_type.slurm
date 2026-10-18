#!/bin/bash
#SBATCH -p short -A stajichlab -c 4 --mem 8gb -t 4:00:00 -J fetch-type
# One SRA run per array task: download to $SCRATCH, type against a species panel, keep only the result row, delete the reads.
# Nothing but the small result TSV lands under /bigdata.
#
# Usage:
#   N=$(( $(wc -l < selected.tsv) - 1 ))
#   sbatch --array=1-${N}%8 fetch_type.slurm <selected.tsv> <downloader.sh> <panel_dir> <panel_prefix> <out_dir>
# e.g. sbatch --array=1-312%8 fetch_type.slurm cryptococcus/selected.tsv cryptococcus/download.sh panels/cryptococcus neoformans out/crypto
# Result file: <out_dir>/<run>.<panel_prefix>.tsv, one per run and panel; datasets/reads/accumulate.py merges them (rerunnable, any order, any number of workers).
# Env: DL_STYLE=outdir_first for cryptococcus/download.sh (default tsv_first, fusarium/download.sh).
# Env: IDIOMORPHS="MAT1-1 MAT1-2" (default) names the idiomorphs; the panel files are <panel_dir>/<panel_prefix>_<NAME>.fasta (e.g. IDIOMORPHS="MATa MATalpha" for Cryptococcus).
# Env: source datasets/reads/coccidioides/pilot/env.sh first (or set CLONE/TOOLS the same way).
set -uo pipefail
: "${SCRATCH:?SCRATCH must be set by Slurm}"
tsv=$1; dl=$2; pdir=$3; pfx=$4; out=$5
run=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n{print $1}' "$tsv")
mkdir -p "$out"
[[ -s $out/$run.$pfx.tsv ]] && { echo "skip $run"; exit 0; }
work=$(mktemp -d "$SCRATCH/fetch_type.$run.XXXXXX"); trap 'rm -rf "$work"' EXIT
# downloader calling conventions: fusarium/download.sh TSV OUTDIR RUN; cryptococcus/download.sh OUTDIR TSV (row = SLURM_ARRAY_TASK_ID, same row as above)
if [[ ${DL_STYLE:-tsv_first} == outdir_first ]]; then bash "$dl" "$work/reads" "$tsv" || { echo "FAIL download $run" >&2; exit 1; }
else bash "$dl" "$tsv" "$work/reads" "$run" || { echo "FAIL download $run" >&2; exit 1; }; fi
# Keep every Nth read pair (about 8M pairs): old SRA FASTQs are position-sorted, so the FIRST N reads can all come from one chromosome.
rc=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="read_count")c=i} NR>1&&$1==r{print $c}' r="$run" "$tsv")
step=$(( (${rc:-0} + 7999999) / 8000000 )); [[ $step -lt 1 ]] && step=1
for m in 1 2; do zcat "$work/reads/$run/${run}_$m.fastq.gz" | awk -v s=$step 'NR%4==1{k=((NR-1)/4)%s==0} k' | gzip > "$work/sub_$m.fastq.gz"; done
printf 'sample\treads\n%s\t%s,%s\n' "$run" "$work/sub_1.fastq.gz" "$work/sub_2.fastq.gz" > "$work/samples.tsv"
idio=(); for n in ${IDIOMORPHS:-MAT1-1 MAT1-2}; do idio+=(--idiomorph "$n=$pdir/${pfx}_$n.fasta"); done
python -m MATPredict reads-type "${idio[@]}" \
  --samples "$work/samples.tsv" --out "$work/$run.$pfx.tsv" || { echo "FAIL type $run" >&2; exit 1; }
mv "$work/$run.$pfx.tsv" "$out/$run.$pfx.tsv"   # result only; the trap removes the reads
