#!/bin/bash
# The standing queries of a sketch store, run after every (re)build: qc, distances (PCoA, clone groups, tree), MAT, and each target
# screen in STORE/targets/screens.tsv (name <TAB> fasta[,fasta...] <TAB> optional extra sketch_query.py screen arguments).
set -uo pipefail
STORE=${1:?STORE}; source "$STORE/env.sh"
PYQ=${PYQ:-/opt/linux/rocky/8.x/x86_64/pkgs/sourmash_branchwater/0.9.13/bin/python}
Q="$PYQ -u $STORE/bin/sketch_query.py --store $STORE --threads ${SLURM_CPUS_PER_TASK:-8}"
echo "== qc";        $Q qc
echo "== mat";       $Q mat --panel "$PANEL_DIR/${LOCAL_LINEAGE:-L1}" --prefix "$PANEL_PREFIX"
if [[ -s $STORE/targets/screens.tsv ]]; then
  while IFS=$'\t' read -r name fastas extra; do
    [[ -z $name || $name == \#* ]] && continue
    echo "== screen $name"
    eval "$Q screen --name $name --targets ${fastas//,/ } $extra"
  done < "$STORE/targets/screens.tsv"
fi
echo "== distances"; $Q distances
echo "done $(date)"
