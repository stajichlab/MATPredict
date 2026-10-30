#!/bin/bash
# Sensitive miniprot of the seven curated Saccharomycotina MATa2 proteins plus P. kudriavzevii MATa2 against each window (HPCC)
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cd "$(dirname "$0")"
for f in data/*.fa; do b=${f%.fa}; miniprot -t4 --gff -P A2 --outs=0.3 --outc=0.2 -k5 -p 0.3 "$f" a2_refs.faa > "$b.a2.gff" 2> "$b.a2.log"; done
