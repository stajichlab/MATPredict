#!/usr/bin/bash -l
#SBATCH --job-name=tremscan
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=48gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/tremello_fail/logs/scan_%A.log
# usage: sbatch --export=ALL,LIST=/path/ids.txt,OUT=/path/outdir scan_all.slurm
source /etc/profile
set -o pipefail
D=/bigdata/stajichlab/jstajich/projects/tremello_fail
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
mkdir -p "$OUT"
run_one() {
  a="$1"
  [[ -s "$OUT/$a.genes.tsv" ]] && return 0
  timeout 1500 /bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python $D/scan/scan_ste3_hd.py "$a" "$OUT" > "$OUT/$a.log" 2>&1 || echo "FAIL $a" >&2
}
export -f run_one
export OUT D
cat "$LIST" | xargs -P "${SLURM_CPUS_PER_TASK:-16}" -I{} bash -c 'run_one {}'
echo done
