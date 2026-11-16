#!/usr/bin/bash -l
#SBATCH --job-name=bd-site
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/bdsite_%A_%a.log
# Depth-controlled read counts at the JEL423/CLFT044 divergent sites. Env: SRC (worktree), BASE (results/2026-11-13_bd_coverage), CHUNK (file, one run per line;
# array task N reads line N when RUNLIST is set), TARGET (depth to downsample to, default 12).
# Steps per run: read-count check against ENA, bwa-mem2 to CLFT044, duplicate removal, downsample to TARGET x (median depth from
# bd_isolate_summary_CLFT044.tsv; no downsampling if the run is at or below TARGET), samtools mpileup at the sites, bd_site_allele_counts.py.
set -uo pipefail
: "${SRC:?}" "${BASE:?}"
TARGET="${TARGET:-12}"
module load sratoolkit/3.2.0 bwa-mem2/2.2.1 samtools/1.22.1
REF="$BASE/ref/CLFT044.fa"; OUT="$BASE/site_counts"; mkdir -p "$OUT"
W="${SCRATCH:?}/bdsite"; mkdir -p "$W"
run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${RUNLIST:?}")
[[ -n "$run" ]] || { echo "no run for task"; exit 1; }
[[ -s "$OUT/$run.counts.tsv.gz" ]] && { echo "$run already done"; exit 0; }
exp=$(awk -F"\t" -v r="$run" '$2==r{print $7}' "$BASE/bd_runs.tsv")
med=$(awk -F"\t" -v r="$run" '$1==r{print $6}' "$BASE/bd_isolate_summary_CLFT044.tsv")
fasterq-dump "$BASE/sra/$run/$run.sra" -e 4 -O "$W" --temp "$W" --split-files -f > "$OUT/$run.fasterq.log" 2>&1 || { echo "FAIL fasterq $run"; exit 1; }
n1=$(( $(wc -l < "$W/${run}_1.fastq") / 4 )); n2=$(( $(wc -l < "$W/${run}_2.fastq") / 4 ))
if [[ "$n1" != "$exp" || "$n2" != "$exp" ]]; then echo "READ COUNT MISMATCH $run mate1=$n1 mate2=$n2 ena=$exp"; exit 1; fi
bwa-mem2 mem -t 4 -R "@RG\tID:$run\tSM:$run" "$REF" "$W/${run}_1.fastq" "$W/${run}_2.fastq" 2> "$OUT/$run.bwa.log" \
  | samtools fixmate -m -@ 2 - - | samtools sort -@ 2 -T "$W/st_$run" - | samtools markdup -r -@ 2 - "$W/$run.bam"
frac=$(awk -v m="$med" -v t="$TARGET" 'BEGIN{ if (m<=t) print 1; else printf "%.4f", t/m }')
if [[ "$frac" == "1" ]]; then cp "$W/$run.bam" "$W/$run.ds.bam"; else samtools view -b -s "7${frac#0}" -@ 2 "$W/$run.bam" > "$W/$run.ds.bam"; fi
samtools index "$W/$run.ds.bam"
echo "$run median_depth=$med downsample_fraction=$frac reads=$n1"
samtools mpileup -l "$BASE/ref/divergent_sites.bed" -f "$REF" -q 20 -Q 20 -B -d 1000 "$W/$run.ds.bam" 2> /dev/null \
  | python3 "$SRC/scripts/bd_site_allele_counts.py" --sites "$BASE/ref/divergent_sites.tsv" --out "$OUT/$run.counts.tsv.gz" | tee "$OUT/$run.summary.txt"
echo "$run done"
