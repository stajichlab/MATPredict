#!/usr/bin/bash -l
#SBATCH --job-name=ffs-aug
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=3:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/ffs_aug3_%A.log
module load hmmer/3.4
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
W=$SCRATCH/ffs; mkdir -p $W
diamond makedb --in /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-13_flank_first_scan/controls_all.faa -d $W/controls --quiet
while read -r g; do
  P=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_gap_lineages2/aug3/$g/predict_results
  python3 /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/flank-first-scan/scripts/flank_first_scan.py --name $g --proteins $P/$g.proteins.fa.gz --gff3 $P/$g.gff3.gz --pfam-hmm /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-13_flank_first_scan/stage1_all/pfam/hmg_homeo.hmm --controls-dmnd $W/controls.dmnd --controls-faa /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-13_flank_first_scan/controls_all.faa --out /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_gap_lineages2/stage1_aug3/per_genome --window 5 --threads 2 || echo FAIL $g
done < /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-14_gap_lineages2/aug3/names.txt
echo done
