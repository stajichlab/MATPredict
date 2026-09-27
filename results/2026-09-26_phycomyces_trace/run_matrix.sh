#!/usr/bin/bash
# Code x input matrix for Phycomyces NRRL_1554. SRC and DB from the same frozen tree.
set -u
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default
export PATH="$E/bin:$PATH"
WT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees
T=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_phycomyces_trace
G=Phycomyces_blakesleeanus_NRRL_1554
CONTIGS=/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/annotate/$G/annotate_results/$G.contigs.fsa
run(){ sha=$1; arm=$2; out=$T/runs/$sha/$arm; mkdir -p $out
  case $arm in
    scaf_prot) args="--genome $T/input/$G.fna --proteins $T/input/$G.faa";;
    scaf)      args="--genome $T/input/$G.fna";;
    contigs)   args="--genome $CONTIGS";;
  esac
  MATPREDICT_DB_ROOT=$WT/run-$sha/db PYTHONPATH=$WT/run-$sha/src timeout 3600 $E/bin/python -m MATPredict detect \
    $args --phylum Mucoromycota --out-dir $out --evidence-diagnostics $out/evidence_diagnostics.jsonl \
    > $out/stdout.log 2> $out/stderr.log; echo "$sha $arm rc=$?"; }
export -f run; export E WT T G CONTIGS
for s in "$@"; do for a in ${ARMS:-scaf_prot scaf contigs}; do echo "$s $a"; done; done | xargs -P 6 -L1 bash -c 'run "$0" "$1"'
