#!/usr/bin/bash -l
#SBATCH --job-name=stringtie
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/stringtie_%A.log
module load stringtie/3.0.3 gffread/0.12.7 hmmer/3.4
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_taphrinomycotina_probe/rnaseq/stringtie
for n in Tdef_A2 Plact_12-1054; do
  stringtie /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_taphrinomycotina_probe/rnaseq/$n/$n.bam -p 4 -o $n.gtf -l ST_$n 2> $n.st.log || echo FAIL $n
  echo "$n transcripts: $(grep -c -P '\ttranscript\t' $n.gtf)"
done
zcat -f /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/GCA_005281805.1_ASM528180v1.fa.gz > $SCRATCH/A2.fa; zcat -f /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/GCF_002105105.1_Prola1.fa.gz > $SCRATCH/L.fa
gffread -w Tdef_A2.tx.fa -g $SCRATCH/A2.fa Tdef_A2.gtf; gffread -w Plact_12-1054.tx.fa -g $SCRATCH/L.fa Plact_12-1054.gtf
echo done
