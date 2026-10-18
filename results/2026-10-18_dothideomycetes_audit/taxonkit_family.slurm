#!/usr/bin/bash -l
#SBATCH -N 1 -c 1 --mem 2gb --time 00:15:00 -J doth-taxfam -p short
# family / order / genus per taxid by taxonkit reformat (never on a login node)
set -euo pipefail
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test/bin:$PATH
cd "$1"
taxonkit lineage taxids_under_147541.txt | taxonkit reformat -i 2 -f "{o}\t{f}\t{g}" -F -r NA > taxid_ranks.tsv
head -3 taxid_ranks.tsv
wc -l taxid_ranks.tsv
