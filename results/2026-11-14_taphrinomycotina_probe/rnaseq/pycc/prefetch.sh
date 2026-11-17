#!/usr/bin/bash -l
module load sratoolkit/3.2.0
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_taphrinomycotina_probe/rnaseq/pycc
mkdir -p sra
while read r n t; do prefetch $r -O sra --max-size 30G > sra/prefetch_$r.log 2>&1; echo "$r rc=$? $(stat -c %s sra/$r/$r.sra 2>/dev/null)"; done < runs.tsv > prefetch_summary.txt 2>&1
echo finished >> prefetch_summary.txt
