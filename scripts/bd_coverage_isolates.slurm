#!/usr/bin/bash
#SBATCH --job-name=bd-cov
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/bdcov_%A_%a.log
# Align Batrachochytrium WGS runs to a reference, then depth windows and SNP allele depths.
# Env: SRC, BASE (results/2026-11-13_bd_coverage), REFNAME (JEL423 or CLFT044; fasta in $BASE/ref), CHUNK (file: one run accession per line).
# Array task N reads $BASE/chunks/chunk_N.txt when CHUNK is unset.
set -uo pipefail
: "${SRC:?}" "${BASE:?}"
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
module load sratoolkit/3.2.0 bwa-mem2/2.2.1 samtools/1.22.1 bcftools/1.22 mosdepth/0.3.12
CH="${CHUNK:-$BASE/chunks/chunk_${SLURM_ARRAY_TASK_ID}.txt}"
REFNAME="${REFNAME:-JEL423}"; REF="$BASE/ref/$REFNAME.fa"
OUT="$BASE/per_run_$REFNAME"; mkdir -p "$OUT"
W="${SCRATCH:?}/bdcov"; mkdir -p "$W"
while read -r run; do
  [[ -z "$run" ]] && continue
  [[ -s "$OUT/$run.snps.tsv.gz" ]] && continue
  [[ -s "$BASE/sra/$run/$run.sra" ]] || { echo "NOT DOWNLOADED $run"; continue; }
  exp=$(awk -F"\t" -v r="$run" '$2==r{print $7}' "$BASE/bd_runs.tsv")
  fasterq-dump "$BASE/sra/$run/$run.sra" -e 8 -O "$W" --temp "$W" --split-files -f > "$OUT/$run.fasterq.log" 2>&1 || { echo "FAIL fasterq $run"; continue; }
  n1=$(( $(wc -l < "$W/${run}_1.fastq") / 4 )); n2=$(( $(wc -l < "$W/${run}_2.fastq") / 4 ))
  if [[ "$n1" != "$exp" || "$n2" != "$exp" ]]; then echo "READ COUNT MISMATCH $run mate1=$n1 mate2=$n2 ena=$exp"; rm -f "$W/${run}"_*.fastq; continue; fi
  bwa-mem2 mem -t 8 -R "@RG\tID:$run\tSM:$run" "$REF" "$W/${run}_1.fastq" "$W/${run}_2.fastq" 2> "$OUT/$run.bwa.log" \
    | samtools fixmate -m -@ 4 - - | samtools sort -@ 4 -T "$W/st_$run" - | samtools markdup -r -@ 4 - "$W/$run.bam"
  samtools index "$W/$run.bam"
  mosdepth -n -t 4 --by 10000 --mapq 20 "$OUT/$run" "$W/$run.bam"
  samtools flagstat "$W/$run.bam" > "$OUT/$run.flagstat.txt"
  bcftools mpileup -f "$REF" -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP -d 1000 -Ou "$W/$run.bam" \
    | bcftools call -mv -V indels -Ou | bcftools view -i 'QUAL>=30' -Ou \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP\t%AD]\n' | gzip > "$OUT/$run.snps.tsv.gz"
  echo "$run reads=$n1 done"
  rm -f "$W/${run}"_*.fastq "$W/$run.bam"*
done < "$CH"
echo done
