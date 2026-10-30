#!/bin/bash
# miniprot of the curated Saccharomycotina MAT/flank proteins (refs.faa) against each locus window. Run on the HPCC.
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cd "$(dirname "$0")"
for f in data/*.fa; do b=${f%.fa}; miniprot -t4 --gff -P BC "$f" refs.faa > "$b.miniprot.gff" 2> "$b.miniprot.log"; done
