#!/bin/bash
# Map the PB4207 models (models_proteins.faa) onto the other four windows with miniprot (splice-aware): identity and intron positions. HPCC.
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cd "$(dirname "$0")"
for f in data/*.fa; do b=${f%.fa}; miniprot -t4 --gff -P L --outs=0.5 "$f" models_proteins.faa > "$b.lift.gff" 2> "$b.lift.log"; done
