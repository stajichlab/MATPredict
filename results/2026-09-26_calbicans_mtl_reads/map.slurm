#!/usr/bin/bash
#SBATCH --job-name=calb-mtl-reads
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/calb_mtl_reads_%A.log
# Do C. albicans assemblies that show one MTL idiomorph hide a heterozygous a/alpha?
# Stream ~2M reads (R1) per isolate from ENA, map to both idiomorphs + a control.
set -uo pipefail
module load minimap2/2.30 samtools/1.22.1
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_calbicans_mtl_reads
W=${SCRATCH:?}/calb; mkdir -p $W
while IFS=$'\t' read -r acc geno run ftp; do
  curl -s --max-time 1200 "https://$ftp" | gzip -dc 2>/dev/null | head -n 8000000 > $W/$run.fq
  echo "$acc $run reads=$(( $(wc -l < $W/$run.fq) / 4 ))"
  minimap2 -t 8 -ax sr $D/mtl_refs.fa $W/$run.fq 2>/dev/null | samtools view -b -q 20 -F 0x904 - | samtools sort -o $W/$run.bam - && samtools index $W/$run.bam
  while IFS=$'\t' read -r c s e name; do
    d=$(samtools depth -a -r "$c:$((s+1))-$e" $W/$run.bam | awk '{t+=$3; n++} END{printf "%.2f", (n? t/n : 0)}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$acc" "$geno" "$run" "$name" "$d"
  done < $D/regions.bed >> $D/depth.tsv
  rm -f $W/$run.fq
done < $D/runs.tsv
echo done
