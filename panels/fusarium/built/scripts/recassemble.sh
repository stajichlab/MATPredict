#!/bin/bash
#SBATCH -p short -A stajichlab -c 16 --mem 48G -t 2:00:00 -J mrx_asm
# usage: sbatch recassemble.sh RUN [ITERATIONS=3]
# 1 recruit read pairs (either mate) by diamond blastx vs the 6-frame-ORF db of all Fusarium panels + MAT proteins
# 2 AAFTF assemble (SPAdes, in the AAFTF apptainer image) on the recruited pairs; 3 pick contigs with blastx hits to the db; map ALL reads (bwa-mem2) to them, recruit mapped pairs; repeat
source /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm/common.sh
run=$1; IT=${2:-3}; T=16
W=$R/work/$run; mkdir -p $W $R/res/asm; cd $W
fetch $run || exit 1
F1=$S/$run/${run}_1.fastq.gz; F2=$S/$run/${run}_2.fastq.gz
for m in 1 2; do
  diamond blastx -d $R/res/ref/recruit -q $S/$run/${run}_$m.fastq.gz --threads $T -e 1e-5 -k 1 --outfmt 6 qseqid sseqid pident length evalue bitscore -o hits_$m.tsv 2>&1 | tail -1
done
cat hits_1.tsv hits_2.tsv | cut -f1 | sort -u > ids_0.txt
wc -l ids_0.txt
cp hits_1.tsv hits_2.tsv $R/res/asm/ 2>/dev/null; for m in 1 2; do cp hits_$m.tsv $R/res/asm/$run.hits_$m.tsv; done; rm -f $R/res/asm/hits_?.tsv
for i in $(seq 0 $IT); do
  seqtk subseq $F1 ids_$i.txt | gzip -1 > rec_${i}_1.fq.gz; seqtk subseq $F2 ids_$i.txt | gzip -1 > rec_${i}_2.fq.gz
  echo "iter $i: $(wc -l < ids_$i.txt) pairs recruited"
  rm -rf asm_$i wd_$i; $AP AAFTF assemble --method spades -1 rec_${i}_1.fq.gz -2 rec_${i}_2.fq.gz -c $T -m 40 -w wd_$i -o asm_$i.fasta > asm_$i.log 2>&1 || { echo "AAFTF assemble failed"; tail asm_$i.log; }
  seqtk seq -L 300 asm_$i.fasta > c_$i.fa
  diamond blastx -d $R/res/ref/recruit -q c_$i.fa --threads $T -e 1e-10 -k 1 --outfmt 6 qseqid sseqid pident length evalue -o chit_$i.tsv 2>&1 | tail -1
  cut -f1 chit_$i.tsv | sort -u > locus_$i.names; echo "iter $i: $(wc -l < locus_$i.names) locus contigs"
  seqtk subseq c_$i.fa locus_$i.names > locus_$i.fa
  [ $i -eq $IT ] && break
  [ -s locus_$i.fa ] || break
  $AP bwa index locus_$i.fa >/dev/null 2>&1
  $AP bwa mem -t $T locus_$i.fa $F1 $F2 2>/dev/null | $AP samtools view -F 4 -q 20 - | cut -f1 | sort -u > mapped_$i.txt
  cat ids_$i.txt mapped_$i.txt | sort -u > ids_$((i+1)).txt
  [ $(wc -l < ids_$((i+1)).txt) -gt 300000 ] && { echo "recruit set exploded (>300k pairs), stopping at iter $i"; break; }
done
cp locus_$i.fa $R/res/asm/$run.locus_contigs.fa; cp c_$i.fa $R/res/asm/$run.contigs_ge300.fa
grep '>' $R/res/asm/$run.locus_contigs.fa
echo done
