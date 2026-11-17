#!/usr/bin/bash -l
#SBATCH --job-name=bd-selfmap
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/bdselfmap_%j.log
# Map one Bd run to its own assembly and write per-100-kb-window depth and read mismatch rate. Read-only analysis job.
# Env: RUN, SRA (path to .sra), FASTA (assembly, .gz ok), EXP (expected read pairs from ENA), OUT (folder), WIN (window bp, default 100000).
# Steps: fasterq-dump, read-count check against EXP (mismatch stops the job), bwa-mem2, duplicate removal, per window: mean depth (samtools bedcov),
# reads with MAPQ >= 20, and edit distance per aligned base (NM tag; primary alignments with MAPQ >= 20).
set -uo pipefail
: "${RUN:?}" "${SRA:?}" "${FASTA:?}" "${EXP:?}" "${OUT:?}"
WIN="${WIN:-100000}"
module load sratoolkit/3.2.0 bwa-mem2/2.2.1 samtools/1.22.1
W="${SCRATCH:?}/selfmap"; mkdir -p "$W" "$OUT"; cd "$W"
zcat -f "$FASTA" > ref.fa; samtools faidx ref.fa; bwa-mem2 index ref.fa > idx.log 2>&1
fasterq-dump "$SRA" --split-files -e 8 -O . -t . > fd.log 2>&1 || { echo "fasterq-dump failed"; exit 1; }
n1=$(( $(wc -l < ${RUN}_1.fastq) / 4 )); n2=$(( $(wc -l < ${RUN}_2.fastq) / 4 ))
[[ "$n1" == "$EXP" && "$n2" == "$EXP" ]] || { echo "READ COUNT MISMATCH mate1=$n1 mate2=$n2 ena=$EXP"; exit 1; }
bwa-mem2 mem -t 8 -R "@RG\tID:$RUN\tSM:$RUN" ref.fa ${RUN}_1.fastq ${RUN}_2.fastq 2> bwa.log | samtools fixmate -m -@ 2 - - | samtools sort -@ 2 -T st - | samtools markdup -r -@ 2 - $RUN.bam
samtools index $RUN.bam
samtools flagstat $RUN.bam > "$OUT/$RUN.flagstat.txt"
awk -v w=$WIN 'BEGIN{OFS="\t"} {n=int(($2-0)/w); } {print}' ref.fa.fai | awk -v w=$WIN -F'\t' 'BEGIN{OFS="\t"} {for(s=0;s<$2;s+=w){e=s+w; if(e>$2)e=$2; print $1,s,e}}' > win.bed
samtools bedcov win.bed $RUN.bam | awk -F'\t' 'BEGIN{print "contig\tstart\tend\tmean_depth"} {print $1"\t"$2"\t"$3"\t"$4/($3-$2)}' > "$OUT/$RUN.depth_windows.tsv"
samtools view -F 0x904 -q 20 $RUN.bam | awk -v w=$WIN -F'\t' 'BEGIN{OFS="\t"} {nm=-1; for(i=12;i<=NF;i++) if($i ~ /^NM:i:/) nm=substr($i,6); if(nm<0) next; cg=$6; al=0; while (match(cg,/^[0-9]+[MIDNSHX=]/)) {n=substr(cg,1,RLENGTH-1)+0; o=substr(cg,RLENGTH,1); cg=substr(cg,RLENGTH+1); if(o=="M"||o=="="||o=="X") al+=n}; k=$3"\t"int($4/w)*w; r[k]++; a[k]+=al; m[k]+=nm} END{print "contig\tstart\treads_q20\taligned_bp\tedit_dist"; for(k in r) print k,r[k],a[k],m[k]}' > "$OUT/$RUN.mismatch_windows.tsv"
echo done
