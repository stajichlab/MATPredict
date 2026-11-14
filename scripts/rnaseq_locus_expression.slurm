#!/usr/bin/bash -l
#SBATCH --job-name=rna-locus
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=4
#SBATCH --mem=24gb
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/rnalocus_%A_%a.log
# Map public RNA-seq to a genome and measure depth over annotated genes and over named loci.
# Env: SRC (worktree), OUT (folder), RUN (SRA accession; .sra already fetched under $OUT/sra/$RUN), ENA_COUNT (read pairs),
#      NAME (label), FASTA, GFF3, LOCI (TSV: name chrom start end, 1-based inclusive; extra regions to report), REPORT (space-separated
#      names to report: gene ids and loci names).
# Steps: fasterq-dump, read-count check against ENA (fail on mismatch), hisat2, samtools bedcov, locus_expression.py.
set -uo pipefail
: "${SRC:?}" "${OUT:?}" "${RUN:?}" "${ENA_COUNT:?}" "${NAME:?}" "${FASTA:?}" "${GFF3:?}" "${LOCI:?}" "${REPORT:?}"
module load sratoolkit/3.2.0 hisat2/2.2.1 samtools/1.22.1
W="${SCRATCH:?}/rna_$NAME"; mkdir -p "$W" "$OUT/$NAME"
cd "$W"
fasterq-dump "$OUT/sra/$RUN/$RUN.sra" --split-files -e ${SLURM_CPUS_PER_TASK:-4} -O "$W" -t "$W" > fd.log 2>&1 || { echo "fasterq-dump failed"; tail -3 fd.log; exit 1; }
n1=$(( $(wc -l < "${RUN}_1.fastq") / 4 )); n2=$(( $(wc -l < "${RUN}_2.fastq") / 4 ))
echo "reads: mate1=$n1 mate2=$n2 ena=$ENA_COUNT"
if [[ "$n1" != "$ENA_COUNT" || "$n2" != "$ENA_COUNT" ]]; then echo "READ COUNT MISMATCH"; exit 1; fi
zcat -f "$FASTA" > genome.fa
hisat2-build -p ${SLURM_CPUS_PER_TASK:-4} genome.fa idx > hb.log 2>&1
hisat2 -p ${SLURM_CPUS_PER_TASK:-4} --dta --max-intronlen 5000 -x idx -1 "${RUN}_1.fastq" -2 "${RUN}_2.fastq" 2> "$OUT/$NAME/hisat2.log" | samtools sort -@ ${SLURM_CPUS_PER_TASK:-4} -o "$NAME.bam" -
samtools index "$NAME.bam"
# BED of all genes (0-based start) plus the named loci
zcat -f "$GFF3" | awk -F'\t' '$3=="gene"{split($9,a,";"); id=a[1]; sub("ID=","",id); print $1"\t"$4-1"\t"$5"\t"id}' > regions.bed
awk -F'\t' '{print $2"\t"$3-1"\t"$4"\t"$1}' "$LOCI" >> regions.bed
samtools bedcov regions.bed "$NAME.bam" > bedcov.tsv
python3 "$SRC/scripts/locus_expression.py" --bedcov bedcov.tsv --report $REPORT --gene-prefix "${GENE_PREFIX:-}" --out "$OUT/$NAME/locus_expression.tsv"
# spliced reads over each named locus (CIGAR with N)
while IFS=$'\t' read -r lname chrom s e; do
  tot=$(samtools view -c "$NAME.bam" "$chrom:$s-$e"); sp=$(samtools view "$NAME.bam" "$chrom:$s-$e" | awk '$6 ~ /N/' | wc -l)
  echo -e "$lname\treads_overlapping=$tot\tspliced_reads=$sp" | tee -a "$OUT/$NAME/spliced.tsv"
done < "$LOCI"
cp "$NAME.bam" "$NAME.bam.bai" "$OUT/$NAME/" 2>/dev/null
echo done
