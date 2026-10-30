#!/bin/bash
#SBATCH -p exfab -A exfab -c 2 --mem 4G -t 00:30:00 -J pf-augustus
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/pf_augustus_%j.log
# Ab initio AUGUSTUS predictions on the four P. fabianii MAT windows (no RNA-seq exists), three yeast species models.
set -euo pipefail
W=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-30_petasospora/work; cd $W
set +eu; source /etc/profile 2>/dev/null; module load augustus/3.5.0; set -eu
for sp in saccharomyces_cerevisiae_S288C candida_guilliermondii pichia_stipitis kluyveromyces_lactis; do
  augustus --species=$sp --genemodel=complete --strand=both --protein=on --codingseq=on --gff3=on windows.fa > aug_$sp.gff3 2> aug_$sp.err || echo "fail $sp"
done
touch ../DONE_AUG
