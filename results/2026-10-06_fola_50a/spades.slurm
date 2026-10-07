#!/bin/bash
#SBATCH -p exfab -A exfab -c 16 --mem 96gb -t 12:00:00 -J fola50a-spades --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_50a/spades_%j.log
# Whole-genome SPAdes of Fola isolate 50a (trimmed reads from N. L.). Work in $SCRATCH, copy back to /bigdata.
set +eu; source /etc/profile; module load spades/4.2.0; set -eu
O=/bigdata/stajichlab/jstajich/projects/fola_reads/50a_spades
spades.py -1 /rhome/nicolel/bigdata/Fola/05_Pre-Processing/50a_R1_trimmed.fastq.gz -2 /rhome/nicolel/bigdata/Fola/05_Pre-Processing/50a_R2_trimmed.fastq.gz --isolate -t 16 -m 90 -o ${SCRATCH:?}/50a_spades
cp -r $SCRATCH/50a_spades/{contigs.fasta,scaffolds.fasta,spades.log,assembly_graph_with_scaffolds.gfa} $O/
echo done
