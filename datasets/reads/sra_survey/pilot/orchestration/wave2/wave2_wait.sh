#!/bin/bash -l
IDS=31226603,31226604,31226605,31226606,31226607,31234601,31234602,31234603,31234604,31234605,31234498,31234719
for i in $(seq 1 480); do n=$(squeue -h -j $IDS | wc -l); [ $n -eq 0 ] && break; sleep 30; done
S=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit; W=/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2
cd $S
{ echo remaining=$n; python3 wave_driver.py gate --wave $W --stage test 2>&1 | tail -1; python3 plan_tally.py --wave $W | cut -c1-200; } > $W/wait_result.txt 2>&1
