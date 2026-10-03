#!/bin/bash
# Direct tblastn (gencode 12) + miniprot of C. auris MTL refs against each genome.
set -euo pipefail
OUT=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_cauris_uncalled
BIN=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
WK=${SCRATCH:?}/cauris_direct; mkdir -p $WK $OUT/direct
for a in $(cat $OUT/genomes.txt); do
  [ -s $OUT/direct/$a.tblastn.tsv ] && continue
  zcat $LIB/$a.fa.gz > $WK/$a.fna
  python3 $OUT/asm_stats.py $WK/$a.fna > $OUT/direct/$a.stats.txt
  $BIN/makeblastdb -in $WK/$a.fna -dbtype nucl -out $WK/$a >/dev/null
  $BIN/tblastn -query $OUT/cauris_refs.faa -db $WK/$a -db_gencode 12 -evalue 1e-5 -max_target_seqs 5 -num_threads 8 \
    -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" > $OUT/direct/$a.tblastn.tsv
  $BIN/miniprot -t8 -T12 --gff -N3 $WK/$a.fna $OUT/cauris_refs.faa 2>/dev/null | grep -E 'mRNA|^##PAF' > $OUT/direct/$a.miniprot.txt || true
  rm -f $WK/$a.fna $WK/$a.n*
done
