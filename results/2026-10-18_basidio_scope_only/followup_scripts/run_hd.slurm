#!/usr/bin/bash
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J fam18-hd -p short
# usage: sbatch --array=1-3%2 run_hd.slurm   (line i of lists/hdscan.txt: acc<TAB>taxid<TAB>fa.gz)
set -uo pipefail
W=/bigdata/stajichlab/jstajich/projects/mat-fam18
cd $W
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test/bin:$PATH
mkdir -p rc/hd
read -r acc taxid fa < <(sed -n "${SLURM_ARRAY_TASK_ID}p" lists/hdscan.txt)
python hd_scan.py $acc $fa rc/hd/$acc.json
