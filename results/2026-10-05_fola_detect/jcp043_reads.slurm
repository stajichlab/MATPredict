#!/bin/bash
#SBATCH -p exfab -A exfab -c 8 --mem 16gb -t 1:00:00 -J fola-jcp043 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/fola-reads-check/results/2026-10-05_fola_detect/jcp043_reads_%j.log
# Read-based MAT typing of JCP043 (SRR28734937): map to the F. oxysporum
# MAT1-1 (AB011379.2) and MAT1-2 (AB011378.1) references, then samtools coverage.
set -e
module load minimap2 samtools 2>/dev/null || true
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cat /bigdata/stajichlab/nicolel/Fola/06_Align/Mating_Types/References/MAT-1.fasta /bigdata/stajichlab/nicolel/Fola/06_Align/Mating_Types/References/MAT-2.fasta > ${SCRATCH:?}/ref.fa
minimap2 -t 8 -ax sr $SCRATCH/ref.fa /bigdata/stajichlab/nicolel/Fola/02_fastq/NCBI_SRAs_retrieved/JCP043_SRR28734937_1.fastq.gz /bigdata/stajichlab/nicolel/Fola/02_fastq/NCBI_SRAs_retrieved/JCP043_SRR28734937_2.fastq.gz 2>/dev/null | samtools view -b -F 4 - | samtools sort -o $SCRATCH/jcp043.bam - && samtools index $SCRATCH/jcp043.bam
samtools coverage $SCRATCH/jcp043.bam > /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/fola-reads-check/results/2026-10-05_fola_detect/JCP043.MAT.coverage.txt
cat /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/fola-reads-check/results/2026-10-05_fola_detect/JCP043.MAT.coverage.txt
