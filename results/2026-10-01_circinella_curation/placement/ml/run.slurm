#!/bin/bash
#SBATCH -p epyc -c 8 --mem 16G --time 8:00:00 -J placeml -o ml/placeml_%j.log
# Placement trees (part 1): 107 tips = label-tree set + 2 out-of-group gains.
# IQ-TREE (ModelFinder, 1000 UFBoot + 1000 SH-aLRT) and RAxML-NG (--all, 200
# bootstraps, ModelFinder's model) on BOTH alignments. Root: Fusarium MAT1-2-1.
source /etc/profile.d/modules.sh
module load iqtree/3.0.1 raxml-ng/2.0.2
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_circinella_curation/placement
O=OUT__MAT1-2-1__5518_3639_MAT_combined
for a in aln_hmm aln_mafft_g50; do
  iqtree3 -s $a.afa -st AA -m MFP -B 1000 -alrt 1000 -T 4 -o $O --prefix ml/iq_$a -redo -quiet &
done
wait
for a in aln_hmm aln_mafft_g50; do
  M=$(grep 'Best-fit model' ml/iq_$a.iqtree | sed 's/.*: //;s/ .*//')
  echo "$a RAxML model from ModelFinder: $M"
  raxml-ng --all --msa $a.afa --model "$M" --bs-trees 200 --threads 4 --seed 1 --outgroup $O --prefix ml/rx_$a --redo > ml/rx_$a.stdout 2>&1 &
done
wait
ls ml/*.treefile ml/*.raxml.support
echo done
