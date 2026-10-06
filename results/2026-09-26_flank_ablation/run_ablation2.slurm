#!/usr/bin/bash
#SBATCH --job-name=flank-ablation2
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/flank_ablation_%A.log
# Flank ablation: 60 Serinales genomes x 5 roster/polish variants, same code.
set -uo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict
E=$R/results/2026-09-26_flank_ablation
SRC=$R/.claude/worktrees/run-882aa01/src
[[ -f $SRC/MATPredict/__init__.py ]] || { echo "no SRC on $(hostname)"; exit 1; }
export PATH=$R/.pixi/envs/default/bin:$PATH PYTHONPATH=$SRC MATPREDICT_CACHE_DIR=$R/.matpredict_cache
W=${SCRATCH:?}/ablation; mkdir -p $W
declare -A DB=( [V4_flanksminiprot_noalpha2]=$E/db/V4_flanksminiprot_noalpha2 [V5_flanksminiprot_alpha2]=$R/.claude/worktrees/run-882aa01/db )
while IFS=$'\t' read -r g t sp; do zcat $R/../../../shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/$g.fa.gz > $W/$g.fna 2>/dev/null || zcat /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/$g.fa.gz > $W/$g.fna; done < $E/genomes.tsv
one() { g=$1; t=$2; v=$3; db=$4; out=$E/runs/$v/$g; [[ -s $out/detection_report.yaml ]] && return 0; mkdir -p $out; t0=$(date +%s)
  if [[ $v == V3m* ]]; then MATPREDICT_DB_ROOT=$db python $E/miniprot_flanks.py detect --genome $W/$g.fna --taxid $t --out-dir $out --evidence-diagnostics $out/ed.jsonl > $out/log 2>&1
  else MATPREDICT_DB_ROOT=$db python -m MATPredict detect --genome $W/$g.fna --taxid $t --out-dir $out --evidence-diagnostics $out/ed.jsonl > $out/log 2>&1; fi
  echo $(( $(date +%s) - t0 )) > $out/wall_seconds; }
export -f one; export E W
for v in "${!DB[@]}"; do while IFS=$'\t' read -r g t sp; do echo "$g $t $v ${DB[$v]}"; done < $E/genomes.tsv; done | xargs -P 16 -L1 bash -c 'one "$0" "$1" "$2" "$3"'
echo "done $(ls $E/runs/*/*/detection_report.yaml | wc -l) / 300"
