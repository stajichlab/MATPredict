#!/bin/bash
#SBATCH -p stajichlab -c 4 --mem 8G -t 2:00:00 -J afum_snpgt -o logs/afum_snpgt.%j.log
# Genotype matrix of the 96 test strains: biallelic SNPs variable among them, one row per site, one column per strain (GT).
set -euo pipefail
module load bcftools
W=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-07_sketch_pcoa_mat
VCF=/bigdata/stajichlab/shared/projects/Population_Genomics/Afumigatus_Global/vcf/Afum_global_v5.All.SNP.combined_selected.vcf.gz
cut -f1 $W/strains.tsv > $W/strain_names.txt
bcftools view --threads 4 -S $W/strain_names.txt -m2 -M2 -v snps -a -c1:nref -Ou $VCF \
  | bcftools view -c1 -Ou | bcftools query -f '[%GT\t]\n' | sed 's/\t$//' > $W/snp_gt.body.tsv
{ bcftools view -S $W/strain_names.txt -h $VCF | bcftools query -l - | paste -sd '\t'; cat $W/snp_gt.body.tsv; } | gzip > $W/snp_gt.tsv.gz
rm $W/snp_gt.body.tsv
echo "sites: $(zcat $W/snp_gt.tsv.gz | tail -n +2 | wc -l)"
