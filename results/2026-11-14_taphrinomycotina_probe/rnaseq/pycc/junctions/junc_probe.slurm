#!/usr/bin/bash -l
#SBATCH --job-name=junc-probe
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_taphrinomycotina_probe/rnaseq/pycc/junctions/junc_%j.log
# One-off probe: where are the spliced-read junctions (and their XS strand) in RHGF01000006.1:50500-53000 for three PYCC libraries.
set -uo pipefail
module load sratoolkit/3.2.0 hisat2/2.2.1 samtools/1.22.1
P=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_taphrinomycotina_probe/rnaseq/pycc
F=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/GCA_005281805.1_ASM528180v1.fa.gz
cd "${SCRATCH:?}"; zcat $F > g.fa; hisat2-build -p 8 g.fa idx > hb.log 2>&1
for R in SRR27065922 SRR27065925 SRR27065919; do
  fasterq-dump $P/sra/$R/$R.sra --split-files -e 8 -O . -t . > fd.log 2>&1 || continue
  hisat2 -p 8 --dta --max-intronlen 5000 -x idx -1 ${R}_1.fastq -2 ${R}_2.fastq 2>/dev/null | samtools sort -@4 -o $R.bam -; samtools index $R.bam
  samtools view $R.bam RHGF01000006.1:50500-53000 | awk -v r=$R '$6 ~ /N/ {x="none"; for(i=12;i<=NF;i++) if($i ~ /^XS:A:/) x=substr($i,6); pos=$4; cg=$6; while (match(cg,/^[0-9]+[MIDNSHX=]/)) {n=substr(cg,1,RLENGTH-1)+0; o=substr(cg,RLENGTH,1); cg=substr(cg,RLENGTH+1); if(o=="N"){print r"\t"pos"\t"pos+n-1"\t"x} if(o=="M"||o=="D"||o=="N"||o=="="||o=="X") pos+=n}}' > $P/junctions/$R.junc.tsv
  rm -f ${R}_*.fastq $R.bam $R.bam.bai
done
