#!/bin/bash
#SBATCH -p batch -J dipo-all -c 2 --mem 4G -t 00:40:00 -o /rhome/jstajich/dipo/log2_%a.txt
P=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/default/bin
mkdir -p ~/dipo/all
s=$(( (SLURM_ARRAY_TASK_ID-1)*10 + 1 )); e=$(( SLURM_ARRAY_TASK_ID*10 ))
for g in $(sed -n "${s},${e}p" ~/dipo/all.txt); do
  [ -s ~/dipo/all/out_$g.tsv ] && continue
  t=$(mktemp -d)
  zcat /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/$g.fa.gz > $t/g.fa || continue
  $P/makeblastdb -in $t/g.fa -dbtype nucl -out $t/db >/dev/null
  $P/tblastn -query ~/dipo/dipo_queries.faa -db $t/db -evalue 1e-3 -num_threads 2 -seg no -outfmt '6 qseqid sseqid pident length qlen bitscore sstart send evalue' > ~/dipo/all/out_$g.tsv
  rm -rf $t
done
