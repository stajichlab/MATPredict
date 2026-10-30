#!/bin/bash -l
S=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit; W=/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2
cd $S
{ python3 wave_driver.py gate --wave $W --stage test 2>&1 | tail -1; python3 wave2_final2.py; python3 wave_driver.py gate --wave $W --stage final 2>&1 | tail -1; python3 wave2_clones.py; } > $W/finish_result.txt 2>&1
