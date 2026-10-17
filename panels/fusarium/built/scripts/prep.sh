#!/bin/bash
#SBATCH -p short -A stajichlab -c 4 --mem 8G -t 1:00:00 -J mrx_prep -o /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm/logs/prep.%j.log
source /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm/common.sh
cd $R/res; mkdir -p ref; cd ref
cat $REPO/db/Ascomycota/Hypocreales/*/proteins.faa > matprot.faa
python3 $R/orfdb.py 40 orf_panels.faa $REPO/panels/fusarium/*_MAT1-[12].fasta $REPO/panels/fusarium_oxysporum_fola/*.fa* 
cat matprot.faa orf_panels.faa > recruit.faa
diamond makedb --in recruit.faa -d recruit -p 4 2>&1 | tail -2
diamond makedb --in matprot.faa -d matprot -p 4 2>&1 | tail -2
# F. andiyazi YN28 assembly (MAT1-1 reference source for the pilot)
for acc in "$@"; do
  curl -sS -o $acc.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA&filename=$acc.zip"
  unzip -o -q $acc.zip -d $acc.dl && cat $acc.dl/ncbi_dataset/data/$acc/*.fna > $acc.fna && rm -rf $acc.zip $acc.dl
  grep -c '>' $acc.fna; blastn -version | head -1; makeblastdb -in $acc.fna -dbtype nucl -out $acc >/dev/null
done
ls -la
