#!/bin/bash -l
# Read-level confirmation of read sets GmG flagged as containing a second organism. For one set (array task = a line of IDS_FILE):
# map its first N pairs to Coccidioides references, keep the pairs where NEITHER mate maps, and classify those with GmG. A real
# contaminant shows as an unmapped fraction well above the clean sets', that GmG names as the contaminant genus; a chance second
# genus does not. Usage:
#   sbatch --array=1-23%6 --export=ALL,IDS=ids.tsv,SAMPLES=samples.tsv,REF=ref_prefix,OUT=dir confirm_flagged.slurm
#   IDS: id<TAB>expected second genus (blank for a negative control);  SAMPLES: GmG sample sheet (id, path, expected), R1 then R2 per id
#   REF: a bowtie2 index prefix (setup: confirm_setup.sh).  PAIRS=8000000.
#SBATCH -p short -c 8 --mem 24gb -t 2:00:00 -J confirm_flagged -o %x.%A_%a.log
set -uo pipefail
module load bowtie2/2.5.3 samtools/1.16 2>/dev/null || true
: "${IDS:?}" "${SAMPLES:?}" "${REF:?}" "${OUT:?}"
GMG=${GMG:-/bigdata/stajichlab/shared/projects/GmG_runs/bin/gmg}
IDX=${GMG_INDEX:-/bigdata/stajichlab/shared/projects/GmG_runs/release/2026.10.5/bfd-2026.10.5-s10000.gmgi}
N=${PAIRS:-8000000}
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDS"); id=$(echo "$line" | cut -f1); exp=$(echo "$line" | cut -f2)
[[ -n $id ]] || { echo "no id on line $SLURM_ARRAY_TASK_ID"; exit 1; }
mkdir -p "$OUT/$id"; [[ -s $OUT/$id/row.tsv ]] && { echo "skip $id"; exit 0; }
w=$(mktemp -d "${SCRATCH:-/tmp}/confirm.$id.XXXXXX"); trap 'rm -rf "$w"' EXIT
r1=$(awk -F'\t' -v i="$id" '$1==i{print $2}' "$SAMPLES" | sed -n 1p); r2=$(awk -F'\t' -v i="$id" '$1==i{print $2}' "$SAMPLES" | sed -n 2p)
{ zcat "$r1" | head -n $((4*N)); } > "$w/r1.fq" 2>/dev/null || true
{ zcat "$r2" | head -n $((4*N)); } > "$w/r2.fq" 2>/dev/null || true
pairs=$(( $(wc -l < "$w/r1.fq") / 4 ))
bowtie2 -p "${SLURM_CPUS_PER_TASK:-8}" --very-sensitive -x "$REF" -1 "$w/r1.fq" -2 "$w/r2.fq" 2> "$OUT/$id/bowtie2.log" \
  | samtools fastq -f 12 -F 0x900 -n -1 "$w/un_1.fq" -2 "$w/un_2.fq" -0 /dev/null -s /dev/null - 2>/dev/null
unm=$(( $(wc -l < "$w/un_1.fq") / 4 ))
aligned=$(grep -o '[0-9.]*% overall alignment rate' "$OUT/$id/bowtie2.log" | cut -d% -f1)
gm=NA
if [[ $unm -ge 2000 ]]; then
  printf 'id\tpath\texpected\n%s\t%s\t%s\n%s\t%s\t%s\n' "$id" "$w/un_1.fq" "$exp" "$id" "$w/un_2.fq" "$exp" > "$w/sheet.tsv"
  "$GMG" run --sample-sheet "$w/sheet.tsv" --index "$IDX" -o "$OUT/$id/gmg" --threads "${SLURM_CPUS_PER_TASK:-8}" > "$OUT/$id/gmg.log" 2>&1 && gm=ok
fi
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$exp" "$pairs" "$unm" "${aligned:-NA}" "$gm" > "$OUT/$id/row.tsv"
echo "done $id: $pairs pairs, $unm unmapped"
