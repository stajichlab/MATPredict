#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict; I=$P/results/sra_reads_typing/wave_inputs
O=$P/results/2026-10-28_taphrinomycotina_mat_candidates; mkdir -p $O/loose
export PATH=$P/.pixi/envs/default/bin:$PATH
python3 $P/results/2026-10-26_v070_submit/taph_mat_candidates.py --db $P/.claude/worktrees/run-b317ece/db --genomes $O/genomes.tsv --genome-dir $I/genomes --out $O/loose --reports-glob "$P/results/2026-10-26_v070_ascomycota/wave_*/small_*/runs/GENOME/detection_report.yaml" --threads 4 --min-identity 0.2 --min-coverage 0.2 > $O/loose/run.log 2>&1
