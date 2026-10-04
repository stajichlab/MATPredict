#!/usr/bin/bash -l
#SBATCH -p exfab -A exfab -c 16 --mem 32gb -t 48:00:00 -J mm_iq_one
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat/tree/logs/%x_%j.log
# One IQ-TREE run (same settings as run_tree.sh) on ALN=hmg|full, or the RAxML-NG
# check (ALN=rx_hmg). Split out of run_tree.sh so the full-length tree fits a 24 h job.
set -eo pipefail
module load iqtree/3.0.1 raxml-ng/2.0.2
set -u
T=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat/tree
A=${ALN:?set ALN}
W=${SCRATCH:?}/iq_$A; mkdir -p $W; cd $W
N=${SLURM_CPUS_PER_TASK:-8}
if [[ $A == rx_hmg ]]; then
  cp $T/hmg.afa .
  M=$(grep "Best-fit model" $T/iq_hmg.iqtree | sed 's/.*: //; s/ .*//')
  raxml-ng --all --msa hmg.afa --model "$M" --bs-trees 200 --threads $N --seed 20261003 --prefix rx_hmg > /dev/null
  cp rx_hmg.* $T/
else
  cp $T/$A.afa .
  # resume from a checkpoint copied back by an earlier run (IQ-TREE reads iq_$A.ckp.gz)
  mkdir -p $T/ckp_$A; cp $T/ckp_$A/iq_$A.* . 2>/dev/null || true
  iqtree3 -s $A.afa --prefix iq_$A -T $N -m MFP -mset LG,WAG,JTT,VT,Q.pfam -B 1000 -alrt 1000 -seed 20261003 > /dev/null &
  pid=$!
  # copy the checkpoint to /bigdata every 30 min: a time limit must not lose the run
  while kill -0 $pid 2>/dev/null; do sleep 1800; cp iq_$A.ckp.gz iq_$A.log iq_$A.model.gz $T/ckp_$A/ 2>/dev/null || true; done
  wait $pid
  cp iq_$A.* $T/
fi
echo done
