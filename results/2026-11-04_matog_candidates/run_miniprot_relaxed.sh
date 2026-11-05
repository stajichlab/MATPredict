#!/bin/bash
# Relaxed miniprot of the two MATog records against each window (run on the HPCC); looks for the short core-gene pieces the default run misses.
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cd "$(dirname "$0")"
for f in win/*.fa; do b=${f%.fa}; miniprot -t4 -k5 -l2 --outs=0.3 --outc=0.2 -N 5 --gff -P MR "$f" refs_og.faa > $b.relaxed.gff 2> $b.relaxed.log; done
