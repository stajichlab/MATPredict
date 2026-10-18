#!/bin/bash
#SBATCH -p short -A stajichlab -c 2 --mem 8G -t 2:00:00 -J mrx_type
# usage: sbatch typejob.sh RUN PREFIX1 [PREFIX2 ...]   (PREFIX = path such that PREFIX_MAT1-1.fasta and PREFIX_MAT1-2.fasta exist)
# downloads RUN, takes every Nth read pair (~1.6M pairs), runs the repo's reads-type code (copied in typer/) once per panel
source /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm/common.sh
run=$1; shift
mkdir -p $R/res/typing
fetch $run || exit 1
subsample $run 1600000
rm -f $S/$run/${run}_?.fastq.gz
# independent truth: diamond blastx of the same subsample vs the MAT proteins (reads per gene)
for m in 1 2; do diamond blastx -d $R/res/ref/matprot -q $S/$run/sub_$m.fastq.gz --threads 2 -e 1e-5 -k 1 --outfmt 6 qseqid sseqid pident -o $S/$run/h$m.tsv >/dev/null 2>&1; done
cat $S/$run/h1.tsv $S/$run/h2.tsv | awk -F"\t" '{split($2,a,"name=");split(a[2],b,"|"); print $1"\t"b[1]"\t"$3}' | sort -u -k1,2 | awk -F"\t" '{n[$2]++; p[$2]+=$3} END{for(g in n) printf "%s:%d(%.0f%%) ", g, n[g], p[g]/n[g]; print ""}' > $R/res/typing/$run.genes.txt
for p in "$@"; do
  python3 $R/typer/type_run.py $run $R/res/typing/$run.$(basename $p).tsv MAT1-1=${p}_MAT1-1.fasta MAT1-2=${p}_MAT1-2.fasta -- $S/$run/sub_1.fastq.gz $S/$run/sub_2.fastq.gz
done
