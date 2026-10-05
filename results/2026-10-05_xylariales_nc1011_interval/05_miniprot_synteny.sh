#!/bin/bash
#SBATCH -p batch -c 32 --mem 64G -t 04:00:00 -J xyl-mpsyn
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/xyl_mpsyn_%j.log
set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/xylariales-interval/results/2026-10-05_xylariales_nc1011_interval/common.sh
set +eu; source /etc/profile 2>/dev/null; module load samtools/1.19.2; set -eu
PY=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python
mkdir -p $D/miniprot; cp $D/synteny/q14.faa $D/miniprot/q14.faa
# HMG-like references: Sordariomycete MAT1-2-1/MAT1-1-3-type proteins from the database + NC1011 KAI0192626.1
for d in $REPO/db/Ascomycota/{Sordariales,Hypocreales,Ophiostomatales,Diaporthales}/*/; do
  [[ -f $d/proteins.faa ]] && awk -v r="$(basename $d)" 'BEGIN{IGNORECASE=1} /^>/{keep=($0 ~ /role=core_MAT/ && $0 ~ /name=(MAT1-2-1|mat a-1|mt a-1|mat1-2-1|MAT1-1-3|matA-3|mat A-3|MAT1-1-3)/); if(keep){n=split($0,a,"name="); split(a[2],b,"|"); print ">DBREF_" r "_" b[1]}; next} keep' $d/proteins.faa
done > $D/miniprot/hmg_refs.faa
awk '/^>/{keep=($1==">KAI0192626.1")} keep' $D/data/ncbi_dataset/data/GCA_022453505.1/protein.faa >> $D/miniprot/hmg_refs.faa
echo "hmg refs: $(grep -c '>' $D/miniprot/hmg_refs.faa)"
cd $D && $PY 05_miniprot_synteny.py $D 14
