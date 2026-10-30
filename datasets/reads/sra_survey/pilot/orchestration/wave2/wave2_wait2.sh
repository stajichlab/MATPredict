#!/bin/bash -l
IDS=31331913,31331914,31331915,31331916,31331917,31331964,31331965,31331966,31331967,31331968
for i in $(seq 1 600); do n=$(squeue -h -j $IDS | wc -l); [ $n -eq 0 ] && break; sleep 30; done
S=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit; W=/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2
cd $S
{ echo remaining=$n; python3 wave_driver.py gate --wave $W --stage test 2>&1 | tail -1; python3 wave2_final.py; python3 wave_driver.py gate --wave $W --stage final 2>&1 | tail -1; } > $W/wait_result2.txt 2>&1
