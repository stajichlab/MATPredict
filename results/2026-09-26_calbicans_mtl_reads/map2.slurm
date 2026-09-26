#!/usr/bin/bash
#SBATCH --job-name=calb-mtl-reads2
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/calb_mtl_reads2_%A.log
# Round 2: WGS runs only, and RANDOM subsampling across the whole R1 file.
# Round 1 took the first 2M reads; some submissions are ordered by genome
# position, so that sampled only the start of chr1 (control 5-8x too deep,
# MTL on chr5 at zero).
set -uo pipefail
module load minimap2/2.30 samtools/1.22.1
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_calbicans_mtl_reads
W=${SCRATCH:?}/calb2; mkdir -p $W
one() { acc=$1; geno=$2; run=$3; ftp=$4; n=$5
  p=$(awk -v n=$n 'BEGIN{p=3000000/n; if(p>1)p=1; printf "%.5f", p}')
  curl -s --max-time 5400 "https://$ftp" | gzip -dc 2>/dev/null | awk -v p=$p -v seed=$RANDOM 'BEGIN{srand(seed)} NR%4==1{keep=(rand()<p)} keep' > $W/$run.fq
  echo "$acc $run sampled_reads=$(( $(wc -l < $W/$run.fq) / 4 )) p=$p"
  minimap2 -t 2 -ax sr $D/mtl_refs.fa $W/$run.fq 2>/dev/null | samtools view -b -q 20 -F 0x904 - | samtools sort -o $W/$run.bam - && samtools index $W/$run.bam
  while IFS=$'\t' read -r c s e name; do
    d=$(samtools depth -a -r "$c:$((s+1))-$e" $W/$run.bam | awk '{t+=$3; k++} END{printf "%.2f", (k? t/k : 0)}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$acc" "$geno" "$run" "$name" "$d"
  done < $D/regions.bed > $W/$run.depth
  rm -f $W/$run.fq
}
export -f one; export D W
xargs -P 4 -L1 bash -c 'one "$0" "$1" "$2" "$3" "$4"' < <(tr '\t' ' ' < $D/runs2.tsv)
cat $W/*.depth > $D/depth2.tsv; echo done
