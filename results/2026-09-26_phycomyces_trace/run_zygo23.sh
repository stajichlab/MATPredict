#!/usr/bin/bash
# All 23 Zygo truth genomes (contigs.fsa, genome-only, --phylum Mucoromycota) at the given commits.
set -u
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default
export PATH="$E/bin:$PATH"
WT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees
T=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_phycomyces_trace
TRUTH=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-23_zygo_bar/zygo_truth.tsv
run(){ sha=$1; org=$2; fna=$3; out=$T/zygo23/$sha/runs/$org
  [[ -s $out/detection_report.yaml ]] && return 0; mkdir -p $out
  MATPREDICT_DB_ROOT=$WT/run-$sha/db PYTHONPATH=$WT/run-$sha/src timeout 3600 $E/bin/python -m MATPredict detect \
    --genome $fna --phylum Mucoromycota --out-dir $out > $out/stdout.log 2> $out/stderr.log || echo "FAIL $sha $org" >&2; }
export -f run; export E WT T
for s in "$@"; do awk -F'\t' -v s=$s '{print s, $1, $6}' $TRUTH | sort -u; done | xargs -P 8 -L1 bash -c 'run "$0" "$1" "$2"'
