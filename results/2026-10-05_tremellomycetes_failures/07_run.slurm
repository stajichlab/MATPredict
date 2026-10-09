#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH -c 16
#SBATCH --mem=32gb
#SBATCH --time=1:00:00
#SBATCH -o /bigdata/stajichlab/jstajich/projects/tremello_fail/logs/inorder_%j.log
D=/bigdata/stajichlab/jstajich/projects/tremello_fail
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
PY=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python
mkdir -p $D/inorder
for k in tri filo cysto; do
  $PY $D/scan/07_inorder_ref_test.py MAKE $D/refs_$k.txt $D/inorder/ref_$k.faa $D/scan_out
  grep -c '>' $D/inorder/ref_$k.faa
  cat $D/list_$k.txt | xargs -P 16 -I{} $PY $D/scan/07_inorder_ref_test.py RUN {} $D/inorder/ref_$k.faa $D/inorder
done
echo done
