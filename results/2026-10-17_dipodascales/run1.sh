#!/bin/bash
#SBATCH -J botrytis-miniprot
#SBATCH -p short
#SBATCH -c 4
#SBATCH -t 1:00:00
#SBATCH --mem 8G
cd /bigdata/stajichlab/jstajich/projects/botrytis-mislabels
mkdir -p mp
G=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
M=/bigdata/stajichlab/jstajich/projects/mat-prof
{ cut -f1 genomes.tsv | sed "s#^#$G/#;s#\$#.fa.gz#"; ls $PWD/ref/*.fna.gz; } > glist.txt
while read f; do
  n=$(basename $f | sed 's/\.fa\.gz$//;s/\.fna\.gz$//')
  zcat $f > mp/$n.fa
  (cd $M && pixi run miniprot -t4 --gff -u -N 3 /bigdata/stajichlab/jstajich/projects/botrytis-mislabels/mp/$n.fa /bigdata/stajichlab/jstajich/projects/botrytis-mislabels/q.faa) > mp/$n.gff 2> mp/$n.err
done < glist.txt
echo done > mp/DONE
