#!/usr/bin/bash -l
module load sratoolkit/3.2.0
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_taphrinomycotina_probe/rnaseq/sra
for r in SRR8132818 SRR5576147; do prefetch $r -O . --max-size 40G > prefetch_$r.log 2>&1; echo "$r rc=$? $(ls -la $r/$r.sra 2>/dev/null | awk '{print $5}')"; done
