#!/bin/bash
#SBATCH -p stajichlab -c 4 --mem 16G -t 2:00:00 -J sketch_calib -o sketch_calib.%j.log
# Calibrate sketch sizes: Af293 genome and a 2M-spot read subsample of one typed A. fumigatus run.
set -euo pipefail
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/sourmash/4.9/bin:/opt/linux/rocky/8.x/x86_64/pkgs/sratoolkit/3.2.0/bin:$PATH
W=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-07_sketch_calibration
mkdir -p $W && cd $W
RUN=SRR13579397
[ -s genome.fna.gz ] || curl -sS -o genome.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz
if [ ! -s ${RUN}_1.fastq.gz ]; then
  fastq-dump -X 2000000 --split-files --gzip $RUN
fi
ls -l genome.fna.gz ${RUN}_*.fastq.gz
for S in 100 1000; do
  sourmash sketch dna -p k=21,k=31,scaled=$S,abund genome.fna.gz -o genome.s$S.sig.zip
  sourmash sketch dna -p k=21,k=31,scaled=$S,abund --merge $RUN ${RUN}_1.fastq.gz ${RUN}_2.fastq.gz -o reads.s$S.sig.zip
  sourmash sig filter -m 2 reads.s$S.sig.zip -o reads.s$S.min2.sig.zip
  sourmash sig filter -m 3 reads.s$S.sig.zip -o reads.s$S.min3.sig.zip
  for k in 21 31; do
    sourmash sig cat -k $k genome.s$S.sig.zip -o genome.s$S.k$k.sig.zip
    sourmash sig cat -k $k reads.s$S.min2.sig.zip -o reads.s$S.min2.k$k.sig.zip
  done
done
for f in *.sig.zip; do
  n=$(sourmash sig fileinfo $f 2>/dev/null | grep -m1 -i 'total hashes' | grep -oE '[0-9]+' | head -1 || true)
  echo -e "SIZE\t$f\t$(stat -c %s $f)\thashes=$n"
done
# containment of genome in reads (completeness of the subsample) and of reads in genome
sourmash sig overlap -k 31 genome.s100.sig.zip reads.s100.min2.sig.zip 2>&1 | tail -12
# MAT1-2 sits in Af293; how many scaled=100 / 1000 hashes fall in a 3.5 kb window? report hash density per Mb instead
echo done
