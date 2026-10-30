#!/bin/bash
# Domain/homology scan of the six-frame segments of each locus window (local; needs pyhmmer, Biopython, Pfam-A.hmm)
cd "$(dirname "$0")"
PY=${PY:-python3}
for f in data/*.fa; do
  b=$(basename "$f" .fa); tag=$(echo "$b" | cut -d_ -f1,2)
  PYTHONPATH=/Users/jstajich/projects/MATPredict/src $PY orf_scan.py "$f" ../2026-10-29_ss_records/refs.faa /Users/jstajich/db/funannotate/Pfam-A.hmm 5000 60000 30 > data/${tag}_orf_scan.tsv 2>&1
done
