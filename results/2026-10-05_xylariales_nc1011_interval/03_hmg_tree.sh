#!/bin/bash
#SBATCH -p batch -c 8 --mem 16G -t 12:00:00 -J xyl-hmgtree
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/xyl_hmgtree_%j.log
# MATA_HMG-family tree: Xylariales HMG proteins among MAT1-2-1/MAT1-1-3 references,
# NCU03481 and fmf-1 (Robinson & Natvig 2019 design, RAxML).
set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/xylariales-interval/results/2026-10-05_xylariales_nc1011_interval/common.sh
set +eu; source /etc/profile 2>/dev/null; module load hmmer/3.4 mafft/7.505 trimal/1.5.1 raxml-ng/1.2.0; set -eu
PY=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python
T=$D/tree; mkdir -p $T/prot; cd $T

# 1. proteomes from NCBI Datasets (labels in headers: label|accession)
tail -n +2 $D/data/tree_genomes.tsv | while IFS=$'\t' read -r acc label grp; do
  [[ -s prot/$label.faa ]] && continue
  if curl -s -m 600 -o prot/$acc.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=PROT_FASTA&filename=x.zip" \
     && unzip -o -q -p prot/$acc.zip "ncbi_dataset/data/$acc/protein.faa" 2>/dev/null | sed "s/^>\([^ ]*\)/>$grp|$label|\1/" > prot/$label.faa && [[ -s prot/$label.faa ]]; then
    echo "ok   $acc $label $(grep -c '>' prot/$label.faa)"
  else
    echo "FAIL $acc $label (no annotated proteins)"; rm -f prot/$label.faa
  fi
  rm -f prot/$acc.zip
done | tee proteome_download.log

# 2. curated Sordariomycete MAT1-2-1 / MAT1-1-3-type references from the database
for d in $REPO/db/Ascomycota/{Sordariales,Hypocreales,Ophiostomatales,Diaporthales}/*/; do
  [[ -f $d/proteins.faa ]] && awk -v r="$(basename $d)" '/^>/{keep=($0 ~ /role=core_MAT/); n=split($0,a,"name="); split(a[2],b,"|"); print ">DBREF|" r "|" b[1]; next} keep' $d/proteins.faa
done > dbref_core.faa
cat prot/*.faa dbref_core.faa > all_proteins.faa
echo "proteins: $(grep -c '>' all_proteins.faa)"

# 3. HMG-box domains
hmmsearch --cpu 8 -E 1e-5 --domtblout all.hmg.domtbl $D/data/hmg_related.hmm all_proteins.faa > /dev/null
$PY $D/extract_hmg.py all_proteins.faa all.hmg.domtbl hmg_domains.faa hmg_domains.tsv

# 4. align, trim, tree
mafft --localpair --maxiterate 1000 --thread 8 --quiet hmg_domains.faa > hmg.aln.faa
trimal -in hmg.aln.faa -out hmg.trim.faa -gt 0.3
raxml-ng --all --msa hmg.trim.faa --model LG+G4 --bs-trees 200 --threads 8 --seed 20261005 --prefix rx --redo
echo done
