#!/usr/bin/env bash
# Submit the raw-read MAT typing of ONE species as separate Slurm jobs that survive a logoff:
#   1. an array, one task per run (stream_type_sra.slurm), at most CONC at a time;
#   2. the same array again after it ends (it skips every finished run, so it only retries failures);
#   3. a summary job after that (summarize_pilot.py) that writes SUMMARY.md and the tables.
# Nothing here needs the submitting shell once the jobs are queued.
#
# Usage: submit_species.sh SLUG RUNS.tsv [--conc 8] [--root DIR] [--target-spots 2000000] [--windows 8] [--dry-run]
#   SLUG      short species name used for the folder and the job names, e.g. pyricularia_oryzae
#   RUNS.tsv  one row per run; columns run_accession, read_count (spots), panel_dir, panel_prefix; optional group, strain,
#             study_accession, library_layout, truth_detect, truth_assembly. panel_dir/panel_prefix point at PANEL_DIR/PREFIX_<IDIOMORPH>.fasta.
# Results (all small; the reads are deleted by each task):
#   ROOT/SLUG/runs.tsv        the input, copied (the record of what was asked)
#   ROOT/SLUG/MANIFEST.txt    when, who, code (git describe), sha256 of the script, run list and panel files, parameters
#   ROOT/SLUG/out/            per run: <run>.<prefix>.tsv (call, breadth, depth, flags) and <run>.run.json (size, seconds, state)
#   ROOT/SLUG/logs/           Slurm logs; ROOT/SLUG/jobs.tsv has the job ids
#   ROOT/SLUG/summary/        SUMMARY.md, summary_by_group.tsv, results_all.tsv, discordant.tsv, failures.tsv
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
slug=${1:?usage: submit_species.sh SLUG RUNS.tsv [options]}; runs=${2:?usage: submit_species.sh SLUG RUNS.tsv [options]}; shift 2
conc=8; root=/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing; target=2000000; win=8; dry=0
while [[ $# -gt 0 ]]; do case $1 in
  --conc) conc=$2; shift 2;; --root) root=$2; shift 2;; --target-spots) target=$2; shift 2;; --windows) win=$2; shift 2;; --dry-run) dry=1; shift;;
  *) echo "unknown option $1" >&2; exit 2;; esac; done
cpus=4; [[ $((conc*cpus)) -le 59 ]] || { echo "conc $conc x $cpus CPUs exceeds the 59 CPU limit (curator ruling: under 60)" >&2; exit 2; }
[[ -s $runs ]] || { echo "cannot read $runs" >&2; exit 2; }
n=$(( $(wc -l < "$runs") - 1 )); [[ $n -ge 1 ]] || { echo "no runs in $runs" >&2; exit 2; }
for c in run_accession read_count panel_dir panel_prefix; do head -1 "$runs" | tr '\t' '\n' | grep -qx "$c" || { echo "RUNS.tsv has no column $c" >&2; exit 2; }; done
d=$root/$slug; py=${PY:-/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python}
sub() { if [[ $dry == 1 ]]; then echo "DRY: sbatch $*" >&2; echo "dry$RANDOM"; else sbatch --parsable "$@"; fi; }
if [[ $dry == 0 ]]; then
  mkdir -p "$d/out" "$d/logs" "$d/summary"; cp "$runs" "$d/runs.tsv"
  { echo "species: $slug"; echo "submitted: $(date -Is) by $USER on $(hostname)"; echo "runs: $n   conc: $conc   target_spots: $target   windows: $win"
    echo "code: $(git -C "$here" describe --tags --always --dirty 2>/dev/null || echo unknown)"
    echo "sha256 script: $(sha256sum "$here/stream_type_sra.slurm" | cut -c1-16)  summarize: $(sha256sum "$here/summarize_pilot.py" | cut -c1-16)  runs: $(sha256sum "$runs" | cut -c1-16)"
    echo "source: https://sra-pub-run-odp.s3.amazonaws.com/sra (NCBI SRA on AWS open data)"
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="panel_dir")a=i;if($i=="panel_prefix")b=i}} NR>1{print $a"\t"$b}' "$runs" | sort -u | while IFS=$'\t' read -r pd pf; do
      for f in "$pd/${pf}"_*.fasta; do echo "panel $(basename "$f"): $(sha256sum "$f" | cut -c1-16) $(grep -v '>' "$f" | tr -d '\n' | wc -c) bp"; done; done; } > "$d/MANIFEST.txt"
fi
env_x="ALL,TARGET_SPOTS=$target,WINDOWS=$win,WORKER_DIR=$here"      # WORKER_DIR: where the Slurm wrapper finds genotype_mat_one_accession.sh (Slurm runs a spooled copy of the wrapper)
common=(-p "${PARTITION:-exfab}" -A "${ACCOUNT:-exfab}" --export="$env_x")
a1=$(sub "${common[@]}" -J "sra-$slug" --array=1-"$n"%"$conc" -o "$d/logs/%x_%A_%a.out" "$here/stream_type_sra.slurm" "$d/runs.tsv" "$d/none" none "$d/out")
a2=$(sub "${common[@]}" -J "sra-$slug-retry" --array=1-"$n"%"$conc" --dependency=afterany:"$a1" -o "$d/logs/%x_%A_%a.out" "$here/stream_type_sra.slurm" "$d/runs.tsv" "$d/none" none "$d/out")
sm=$(sub -p "${PARTITION:-exfab}" -A "${ACCOUNT:-exfab}" -c 1 --mem 2gb -t 0:15:00 -J "sra-$slug-summary" --dependency=afterany:"$a2" -o "$d/logs/%x_%j.out" \
     --wrap "$py $here/summarize_pilot.py --runs $d/runs.tsv --results $d/out --out $d/summary")
[[ $dry == 1 ]] || printf 'array\t%s\nretry\t%s\nsummary\t%s\n' "$a1" "$a2" "$sm" > "$d/jobs.tsv"
echo "$slug: $n runs, $conc at a time. jobs: array $a1, retry $a2, summary $sm. Results: $d (summary in $d/summary/SUMMARY.md)"
echo "check: sacct -j $a1,$a2,$sm -X -o JobID,JobName,State,Elapsed ; ls $d/out/*.tsv | wc -l"
