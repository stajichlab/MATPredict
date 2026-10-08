#!/bin/bash
#SBATCH -p short -A stajichlab -c 4 --mem 8gb -t 2:00:00 -J afum-blast --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/blast/slurm_%j.log
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
mkdir -p ${SCRATCH:?}/bt
cat /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/genomes.txt | xargs -P4 -I{} /usr/bin/python3.12 /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/afum-reads-type/scripts/assembly_idiomorph_blast.py /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/afum-reads-type/panels/aspergillus_fumigatus/specific_regions.fasta {} $SCRATCH/bt > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/blast/assembly_truth.rows.tsv
(echo -e "strain\tMAT1-1_cov\tMAT1-1_pid\tMAT1-1_state\tMAT1-2_cov\tMAT1-2_pid\tMAT1-2_state"; sort /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/blast/assembly_truth.rows.tsv) > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/blast/assembly_truth.tsv
echo done
