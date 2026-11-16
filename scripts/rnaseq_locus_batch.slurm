#!/usr/bin/bash -l
#SBATCH --job-name=rna-batch
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/rnabatch_%A_%a.log
# Several RNA-seq runs of one strain mapped to one genome; one table row per run. Read-only discovery job.
# Env: SRC (worktree), OUT (folder; .sra files under $OUT/sra/<RUN>/<RUN>.sra), RUNS (TSV: run, ena_read_pairs, label), FASTA, GFF3,
#      LOCI (TSV name, contig, start, end; 1-based inclusive), REPORT (names to report), GENE_PREFIX, SPLIT_REGIONS (space list name:contig:start:end for splice-strand counts),
#      CHUNK (runs per array task; the task takes lines SLURM_ARRAY_TASK_ID*CHUNK+1 ..).
# Per run: fasterq-dump, read-count check against ENA (a mismatch skips the run and is logged), hisat2 --dta, depth percentile of the loci among
# annotated genes (locus_expression.py), strandedness (infer_strandedness.py), spliced reads by splice-motif strand (XS tag) in each split region.
set -uo pipefail
: "${SRC:?}" "${OUT:?}" "${RUNS:?}" "${FASTA:?}" "${GFF3:?}" "${LOCI:?}" "${REPORT:?}" "${CHUNK:?}"
module load sratoolkit/3.2.0 hisat2/2.2.1 samtools/1.22.1
T=${SLURM_ARRAY_TASK_ID:-0}
W="${SCRATCH:?}/rnab_$T"; mkdir -p "$W" "$OUT/rows"; cd "$W"
zcat -f "$FASTA" > genome.fa
hisat2-build -p ${SLURM_CPUS_PER_TASK:-8} genome.fa idx > hb.log 2>&1
zcat -f "$GFF3" | awk -F'\t' '$3=="gene"{split($9,a,";"); id=a[1]; sub("ID=","",id); print $1"\t"$4-1"\t"$5"\t"id}' > genes.bed
awk -F'\t' '{print $2"\t"$3-1"\t"$4"\t"$1}' "$LOCI" > loci.bed
cat genes.bed loci.bed > regions.bed
sed -n "$((T*CHUNK+1)),$(((T+1)*CHUNK))p" "$RUNS" | while IFS=$'\t' read -r RUN ENA LABEL; do
  [[ -z "$RUN" ]] && continue
  fasterq-dump "$OUT/sra/$RUN/$RUN.sra" --split-files -e ${SLURM_CPUS_PER_TASK:-8} -O "$W" -t "$W" > fd_$RUN.log 2>&1 || { echo "$RUN fasterq-dump failed" >> "$OUT/rows/problems.txt"; continue; }
  n1=$(( $(wc -l < "${RUN}_1.fastq") / 4 )); n2=$(( $(wc -l < "${RUN}_2.fastq") / 4 ))
  if [[ "$n1" != "$ENA" || "$n2" != "$ENA" ]]; then echo "$RUN READ COUNT MISMATCH mate1=$n1 mate2=$n2 ena=$ENA" >> "$OUT/rows/problems.txt"; rm -f ${RUN}_*.fastq; continue; fi
  hisat2 -p ${SLURM_CPUS_PER_TASK:-8} --dta --max-intronlen 5000 -x idx -1 "${RUN}_1.fastq" -2 "${RUN}_2.fastq" 2> hs_$RUN.log | samtools sort -@ 4 -o $RUN.bam -
  samtools index $RUN.bam
  rate=$(grep "overall alignment rate" hs_$RUN.log | awk '{print $1}')
  samtools bedcov regions.bed $RUN.bam > bedcov_$RUN.tsv
  python3 "$SRC/scripts/locus_expression.py" --bedcov bedcov_$RUN.tsv --report $REPORT --gene-prefix "${GENE_PREFIX:-}" --out le_$RUN.tsv > /dev/null
  st=$(python3 "$SRC/scripts/infer_strandedness.py" --bam $RUN.bam --gff3 "$GFF3" | head -1)
  line="$RUN\t$LABEL\t$ENA\t$rate"
  for tag in $REPORT; do d=$(awk -F'\t' -v n=$tag '$1==n{print $3"\t"$4}' le_$RUN.tsv); line="$line\t$d"; done
  for reg in $SPLIT_REGIONS; do
    IFS=: read -r rn rc rs re <<< "$reg"
    xs=$(samtools view $RUN.bam $rc:$rs-$re | awk '$6 ~ /N/ {x="none"; for(i=12;i<=NF;i++) if($i ~ /^XS:A:/) x=substr($i,6); c[x]++} END{printf "%d\t%d", c["+"], c["-"]}')
    line="$line\t$xs"
  done
  printf "$line\t$st\n" > "$OUT/rows/$RUN.row.tsv"
  rm -f ${RUN}_*.fastq $RUN.bam $RUN.bam.bai
done
echo done
