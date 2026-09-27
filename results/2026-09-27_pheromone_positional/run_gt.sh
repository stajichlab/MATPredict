#!/bin/bash
# ground-truth genomes: ASMID and same-genus record prefixes to exclude from Hx
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_pheromone_positional
while read asm ex; do $E scan_genome.py $asm out_gt $ex; done <<LIST
GCA_016772295.1_ASM1677229v1 5346_
GCF_000143185.2_Schco3 5334_
GCF_000091045.1_ASM9104v1 40410_,178876_,5207_,37769_
GCF_000328475.2_Umaydis521_2.0 5270_
GCA_921037615.3_Hybrid_genome_assembly_and_annotation_02 5286_
GCA_000988875.2_ASM98887v2 5286_
LIST
