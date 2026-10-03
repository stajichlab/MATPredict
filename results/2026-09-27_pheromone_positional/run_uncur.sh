#!/bin/bash
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_pheromone_positional
cut -f1 uncurated_sample.tsv | while read a; do [ -s out_uncur/$a.loci.tsv ] || $E scan_genome.py $a out_uncur 2>>out_uncur.err; done
echo done
