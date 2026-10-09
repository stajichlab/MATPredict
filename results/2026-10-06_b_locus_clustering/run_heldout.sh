#!/bin/bash
#SBATCH -p short -c 8 --mem 16G -t 1:00:00 -J bl_hx
# run from the worktree copy on HPCC scratch; reads genomes only
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
export SCRATCH=${SCRATCH:-/tmp}
PY=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python
O=$PWD/hx_out
while read asm excl; do $PY heldout_hx.py $asm $excl $O; done <<EOL
GCF_000143185.2_Schco3 5334_
GCA_016772295.1_ASM1677229v1 5346_
GCF_000271585.1_Trametes_versicolor_v1.0 none
GCF_000320585.1_Heterobasidion_irregulare_v2.0 none
GCA_001683735.1_ASM168373v1 none
GCA_984573805.1_gfRusNobi1.hap1.1 none
EOL
