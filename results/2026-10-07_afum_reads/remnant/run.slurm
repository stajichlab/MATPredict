#!/bin/bash
#SBATCH -p short -A stajichlab -c 8 --mem 24gb -t 2:00:00 -J afum-remnant --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/remnant/slurm_%j.log
# Remnant of the A1163 MAT1-2-1: BLAST + minimap2 liftover onto the A1163 assembly, for all 304 assemblies.
set +eu; source /etc/profile; module load minimap2/2.30; set -eu
export PATH=$PATH:/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
G=/bigdata/stajichlab/shared/projects/Afumigatus_pangenome/scaffolded/genomes
mkdir -p ${SCRATCH:?}/rt
cat /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/genomes.txt | xargs -P4 -I{} /usr/bin/python3.12 /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/afum-reads-type/scripts/afum_remnant.py $G/A1163.sorted.fasta {} $SCRATCH/rt/$(basename {} .sorted.fasta) > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/remnant/rows.raw.tsv
grep -v '^STRIP' /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/remnant/rows.raw.tsv | sort > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/remnant/body.tsv
(echo -e "strain\tremnant_cov\tremnant_pid\tremnant_hit\tmat111_cov\tmat111_hit\tremnant_mat111_gap\tlift_cov\tlift_interval\tsame_location\tmat111_lift_cov"; cat /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/remnant/body.tsv) > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/remnant/remnant_summary.tsv
grep '^STRIP' /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/remnant/rows.raw.tsv | cut -f2- | sort > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/remnant/locus_strip.tsv
echo done
