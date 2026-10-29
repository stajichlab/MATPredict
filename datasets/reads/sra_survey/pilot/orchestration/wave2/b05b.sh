#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict; R=$P/results/sra_reads_typing; W=$R/waves/wave2/plan/botrytis_cinerea/panel
export PATH=$P/.pixi/envs/default/bin:$PATH
T=$(mktemp -d /tmp/b05.XXXX); gunzip -c $R/wave_inputs/genomes/GCF_000143535.2_ASM14353v4.fa.gz > $T/b05.fa
for f in MAT1-1 MAT1-2; do
  echo "== L1/botrytis_$f.fasta vs B05.10 (GCF_000143535.2)"
  grep ">" $W/L1/botrytis_$f.fasta | cut -c1-110
  blastn -query $W/L1/botrytis_$f.fasta -subject $T/b05.fa -evalue 1e-5 -outfmt "6 qseqid qlen length pident qstart qend sseqid sstart send evalue" -max_hsps 6 2>/dev/null | awk -F'\t' '{printf "%s qlen=%s aln=%s id=%.1f q=%s-%s -> %s:%s-%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9}' | sort -k1,1 | head -24
done
