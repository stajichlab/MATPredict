#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict; R=$P/results/sra_reads_typing; S=$P/results/2026-10-26_v070_submit
export PATH=$P/.pixi/envs/default/bin:$PATH
which minimap2 samtools bwa 2>&1 | head -3
O=$R/waves/wave2/botrytis/b05map; mkdir -p $O/out; cd $O
FASTQ_OUT=$O/SRR29713116.fq SKIP_TYPE=1 SCRATCH=/tmp bash $S/genotype_mat_one_accession.sh SRR29713116 34705289 $R/waves/wave2/plan/botrytis_cinerea/panel botrytis $O/out b05 PAIRED 2>&1 | tail -2
ls -la $O | head
