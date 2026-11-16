#!/usr/bin/bash -l
#SBATCH --job-name=plact-hd
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/placthd_%j.log
# Exon-level check of the homeodomain-like region in Protomyces lactucae-debilis 12-1054 with the reverse-stranded library SRR5576147. Read-only.
# Env: SRC (worktree), SRA, FASTA, EXP (ENA read pairs), OUT, HMM (HMG/homeodomain Pfam HMM file), CONTIG, REG_START, REG_END (the region), PAD (default 3000).
# Steps: fasterq-dump, count check, hisat2 --rna-strandness RF --dta, samtools view of the padded region: spliced-read junctions with XS strand and
# counts per junction; strand-resolved depth in 100-bp bins (read 1 reverse = gene strand, dUTP library; checked by infer_strandedness.py in an earlier job);
# StringTie --rf on the whole BAM; the transcripts overlapping the region are extracted, translated in six frames and searched with hmmsearch (E <= 1).
set -uo pipefail
: "${SRC:?}" "${SRA:?}" "${FASTA:?}" "${EXP:?}" "${OUT:?}" "${HMM:?}" "${CONTIG:?}" "${REG_START:?}" "${REG_END:?}"
PAD="${PAD:-3000}"; RUN=$(basename "$SRA" .sra)
module load sratoolkit/3.2.0 hisat2/2.2.1 samtools/1.22.1 stringtie/3.0.3 gffread/0.12.7 hmmer/3.4
W="${SCRATCH:?}/plact_hd"; mkdir -p "$W" "$OUT"; cd "$W"
zcat -f "$FASTA" > g.fa; samtools faidx g.fa; hisat2-build -p 8 g.fa idx > hb.log 2>&1
fasterq-dump "$SRA" --split-files -e 8 -O . -t . > fd.log 2>&1 || { echo "fasterq-dump failed"; exit 1; }
n1=$(( $(wc -l < ${RUN}_1.fastq) / 4 )); n2=$(( $(wc -l < ${RUN}_2.fastq) / 4 ))
[[ "$n1" == "$EXP" && "$n2" == "$EXP" ]] || { echo "READ COUNT MISMATCH mate1=$n1 mate2=$n2 ena=$EXP"; exit 1; }
hisat2 -p 8 --dta --rna-strandness RF --max-intronlen 5000 -x idx -1 ${RUN}_1.fastq -2 ${RUN}_2.fastq 2> hs.log | samtools sort -@ 4 -o $RUN.bam -
samtools index $RUN.bam; grep "overall alignment rate" hs.log > "$OUT/alignment_rate.txt"
S=$((REG_START-PAD)); E=$((REG_END+PAD)); R="$CONTIG:$S-$E"
# junctions (every spliced read segment) with the hisat2 XS tag
samtools view $RUN.bam "$R" | awk '$6 ~ /N/ {x="none"; for(i=12;i<=NF;i++) if($i ~ /^XS:A:/) x=substr($i,6); pos=$4; cg=$6; while (match(cg,/^[0-9]+[MIDNSHX=]/)) {n=substr(cg,1,RLENGTH-1)+0; o=substr(cg,RLENGTH,1); cg=substr(cg,RLENGTH+1); if(o=="N") print pos"\t"pos+n-1"\t"x; if(o=="M"||o=="D"||o=="N"||o=="="||o=="X") pos+=n}}' | sort | uniq -c | awk 'BEGIN{OFS="\t";print "reads\tintron_start\tintron_end\tXS"} {print $1,$2,$3,$4}' > "$OUT/junctions.tsv"
# strand-resolved depth: flag 0x10 and read 1 -> read 1 on the minus strand of the genome means gene on the plus strand for a dUTP library
# mate 1 reverse (0x10) or mate 2 forward = gene on the plus strand (dUTP library); the reverse for the minus strand
samtools view -b -f 0x40 -f 0x10 $RUN.bam "$R" > p1.bam; samtools view -b -f 0x80 -F 0x10 $RUN.bam "$R" > p2.bam
samtools view -b -f 0x40 -F 0x10 $RUN.bam "$R" > m1.bam; samtools view -b -f 0x80 -f 0x10 $RUN.bam "$R" > m2.bam
samtools merge -f plus.bam p1.bam p2.bam; samtools merge -f minus.bam m1.bam m2.bam; samtools index plus.bam; samtools index minus.bam
awk -v c=$CONTIG -v s=$S -v e=$E 'BEGIN{OFS="\t"; for(a=s;a<e;a+=100) print c,a,a+100}' > bins.bed
samtools bedcov bins.bed plus.bam | awk -F'\t' 'BEGIN{print "contig\tstart\tend\tplus_depth"}{print $1"\t"$2"\t"$3"\t"$4/100}' > plus.tsv
samtools bedcov bins.bed minus.bam | awk -F'\t' '{print $4/100}' > minus.col
paste plus.tsv <(echo minus_depth; cat minus.col) > "$OUT/strand_depth_100bp.tsv"
# StringTie, reverse-stranded
stringtie $RUN.bam --rf -p 8 -o st.gtf -l PL > st.log 2>&1
awk -F'\t' -v c=$CONTIG -v s=$S -v e=$E '$1==c && $4<=e && $5>=s' st.gtf > "$OUT/region_stringtie.gtf"
gffread -w "$OUT/region_tx.fa" -g g.fa "$OUT/region_stringtie.gtf" 2> gffread.log
python3 "$SRC/scripts/tx_six_frame.py" --fasta "$OUT/region_tx.fa" --out "$OUT/region_tx_6frame.faa"
hmmsearch -E 1 --domtblout "$OUT/region_tx_6frame.dom" -o /dev/null "$HMM" "$OUT/region_tx_6frame.faa"
echo done
