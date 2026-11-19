set -e
R=/bigdata/stajichlab/jstajich/projects/MATPredict
B=$R/results/2026-11-06_seg_followups
W=$R/.claude/worktrees/seg-switch; F=$R/.claude/worktrees/run-50bc950
cut -f1,2 $B/bas_panel.tsv > $B/lists/bas.tsv
cut -f1,2 $B/muc_panel.tsv > $B/lists/muc.tsv
cut -f1,2 $R/results/2026-11-06_seg_call_test/panel300.tsv | shuf -n 60 --random-source=<(yes 7) > $B/lists/asco60.tsv
wc -l $B/lists/*.tsv
cat > $B/run_muc.slurm <<'S'
#!/bin/bash
#SBATCH -p stajichlab --constraint=milan --time=3:00:00 --cpus-per-task=16 --mem=48gb -J segfu-muc
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-06_seg_followups/muc_%j.log
set -u
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict; W=$REPO/.claude/worktrees/seg-switch; F=$REPO/.claude/worktrees/run-50bc950
B=$REPO/results/2026-11-06_seg_followups; OUT=$B/muc; mkdir -p $OUT
export PYTHONPATH=$W/src MATPREDICT_DB_ROOT=$F/db MATPREDICT_CACHE_DIR=$REPO/.matpredict_cache MATPREDICT_TBLASTN_SEG=yes
PY=$REPO/.pixi/envs/default/bin/python; export PATH=$REPO/.pixi/envs/default/bin:$PATH
S=${SCRATCH:?}/segfu_muc; mkdir -p $S
run_one(){ id="$1"; g="$2"; out="$OUT/$id"; [[ -s "$out/detection_report.yaml" ]] && return 0
  mkdir -p "$out"; fna="$S/$id.fna"; if [[ "$g" == *.gz ]]; then zcat "$g" > "$fna"; else cp "$g" "$fna"; fi
  t0=$(date +%s)
  timeout 3600 "$PY" -m MATPredict detect --genome "$fna" --phylum Mucoromycota --preset v0.8.1 --out-dir "$out" --emit-cds-fasta > "$out/stdout.log" 2> "$out/stderr.log"
  echo "$(( $(date +%s) - t0 ))" > "$out/wall_seconds"; rm -f "$fna"; }
export -f run_one; export OUT S PY
awk -F'\t' '{print $1" "$2}' $B/lists/muc.tsv | xargs -P 16 -L1 bash -c 'run_one "$0" "$1"'
echo "done $(find $OUT -name detection_report.yaml -size +0 | wc -l)"
S
cd $R
sub(){ # name list outdir nodelist-or-none seg dep
 extra=""; [ -n "$4" ] && extra="--nodelist=$4"; [ -n "$6" ] && extra="$extra --dependency=afterok:$6"
 sbatch --parsable --partition=stajichlab --constraint=milan $extra --time=3:00:00 --cpus-per-task=16 --mem=32gb --job-name=segfu-$1 --output=$R/logs/segfu-$1_%j.log \
  --export=ALL,REPO=$R,SRC=$W/src,MATPREDICT_DB_ROOT=$F/db,CLADE=segfu-$1,SAMPLE_LIST=$B/lists/$2.tsv,OUT_DIR=$B/$3,GENOME_TIMEOUT=3600,DETECT_ARGS="--preset v0.8.1",MATPREDICT_TBLASTN_SEG=$5 $F/scripts/run_clade_panel.slurm; }
echo bas $(sub bas bas bas "" yes "")
echo muc $(sbatch --parsable $B/run_muc.slurm)
c1=$(sub t_no asco60 t_no r12 no ""); echo timing_no $c1
echo timing_yes $(sub t_yes asco60 t_yes r12 yes $c1)
