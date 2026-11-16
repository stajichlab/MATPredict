#!/usr/bin/bash -l
#SBATCH --job-name=st-pycc
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/st_pycc_%j.log
# StringTie on the PYCC 5710 24H1 BAM (unstranded; hisat2 --dta gives XS from splice motifs), then translate transcripts at the A2 matPi region.
module load stringtie/3.0.3 gffread/0.12.7 hmmer/3.4
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_taphrinomycotina_probe
H=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-13_flank_first_scan/stage1_all/pfam/hmg_homeo.hmm
cd $D/rnaseq/stringtie_pycc
stringtie $D/rnaseq/Tdef_PYCC5710_24H1/Tdef_PYCC5710_24H1.bam -p 2 -o pycc.gtf -l ST_PYCC 2> pycc.st.log
zcat -f /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/GCA_005281805.1_ASM528180v1.fa.gz > $SCRATCH/A2.fa
gffread -w pycc.tx.fa -g $SCRATCH/A2.fa pycc.gtf
awk -F'\t' '$1=="RHGF01000006.1" && $4<53500 && $5>50500' pycc.gtf > locus.gtf
gffread -y locus.faa -g $SCRATCH/A2.fa locus.gtf
hmmsearch --cut_ga --domtblout locus.domtbl -o /dev/null $H locus.faa
echo done
