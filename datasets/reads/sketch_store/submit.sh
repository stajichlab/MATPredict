#!/bin/bash
# Queue a sketch-store build on the UCR HPCC as a self-feeding planner: each round queues the samples without a sketch that fit under
# the account limits, split over preempt and stajichlab, then queues itself again after them; when nothing is left (or only samples
# that failed MAX_TRIES times), it queues the queries (run_queries.sh).
#
#   submit.sh STORE [AFTER_JOB_IDS]       starts the first round after AFTER_JOB_IDS (colon-separated) finish
#   PLAN=1 submit.sh STORE                (internal) one planner round
#
# STORE must hold manifest.tsv, env.sh (SPECIES, PANEL_DIR, PANEL_PREFIX, LOCAL_INPUT, GENOME_KB, COVERAGE), bin/ (sketch_one.sh,
# sketch_query.py, run_queries.sh, submit.sh, genotype_mat_one_accession.sh, screen_reads.py) and targets/ for the screens.
# Cluster limits (2026-11): array indices below MaxArraySize 2500; 5000 queued array tasks per account (other work shares it). So each
# array is at most CHUNK tasks reading their manifest row numbers from a list file, and a round takes only the room left (minus MARGIN).
# Every attempt is logged in STORE/lists/attempts.tsv; jobs in STORE/jobs.tsv.
set -euo pipefail
STORE=$(readlink -f "${1:?STORE}"); AFTER=${2:-}
source "$STORE/env.sh"
CHUNK=${CHUNK:-2000}; MARGIN=${MARGIN:-300}; MAX_TRIES=${MAX_TRIES:-3}; MAX_ROUNDS=${MAX_ROUNDS:-30}
mkdir -p "$STORE/logs" "$STORE/lists"; touch "$STORE/lists/attempts.tsv" "$STORE/jobs.tsv"
plan_job() {   # plan_job DEPENDENCY: queue the next planner round
  sbatch --parsable ${1:+--dependency=$1} -p stajichlab -A stajichlab -c 1 --mem 2G -t 0:30:00 -J "sk-$SPECIES-plan" \
    -o "$STORE/logs/plan.%j.log" --wrap "PLAN=1 bash $STORE/bin/submit.sh $STORE"
}

if [[ ${PLAN:-0} != 1 ]]; then
  p=$(plan_job "${AFTER:+afterany:$AFTER}")
  echo "planner $p${AFTER:+ (after $AFTER)}"; echo -e "$(date +%F_%T)\tstart\t$p\t$AFTER" >> "$STORE/jobs.tsv"; exit 0
fi

round=$(( $(grep -c $'\tround' "$STORE/jobs.tsv" || true) + 1 ))
# cache index for other tools: every sample with a sketch, its run and the sketch parameters (check here before downloading a run)
awk -F'\t' -v d="$STORE/sigs/" -v p="${PARAMS:-k=21,k=31,scaled=100,abund};min_abund=2;subsample=${COVERAGE:-30}x_genome;sketcher=sourmash_branchwater_singlesketch" \
  'BEGIN{OFS="\t"; print "run_accession","sample_id","source","params","path"} NR>1{f=d $1 ".sig.zip"; if ((getline x < f) > 0) print $3, $1, $2, p, f; close(f)}' \
  "$STORE/manifest.tsv" > "$STORE/sketched.tsv.tmp" && mv "$STORE/sketched.tsv.tmp" "$STORE/sketched.tsv"
# samples without a sketch, minus those already tried MAX_TRIES times
tail -n +2 "$STORE/manifest.tsv" | cut -f1 | awk -v d="$STORE/sigs/" '{print NR "\t" d $1 ".sig.zip"}' \
  | while IFS=$'\t' read -r i f; do [[ -s $f ]] || echo "$i"; done > "$STORE/lists/missing.now"
awk -v m="$MAX_TRIES" 'FILENAME==ARGV[1]{n[$1]++; next} n[$1] < m' "$STORE/lists/attempts.tsv" "$STORE/lists/missing.now" > "$STORE/lists/todo.now"
ntodo=$(wc -l < "$STORE/lists/todo.now"); nmiss=$(wc -l < "$STORE/lists/missing.now")
if [[ $ntodo -eq 0 || $round -gt $MAX_ROUNDS ]]; then
  j=$(sbatch --parsable -p stajichlab -A stajichlab -c "${QUERY_CPUS:-16}" --mem "${QUERY_MEM:-128G}" -t 12:00:00 -J "sk-$SPECIES-query" \
      -o "$STORE/logs/query.%j.log" --wrap "source $STORE/env.sh; bash $STORE/bin/run_queries.sh $STORE")
  echo "sketching finished after $((round - 1)) rounds: $nmiss without a sketch (lists/missing.now); queries $j"
  echo -e "$(date +%F_%T)\tqueries\t$j\t$nmiss missing" >> "$STORE/jobs.tsv"; exit 0
fi
room() { echo $(( 5000 - MARGIN - $(squeue -h -r -u "$USER" -A "$1" 2>/dev/null | wc -l) )); }
rp=$(room preempt); rl=$(room stajichlab); [[ $rp -lt 0 ]] && rp=0; [[ $rl -lt 0 ]] && rl=0
[[ $rp -gt $CHUNK ]] && rp=$CHUNK; [[ $rl -gt $CHUNK ]] && rl=$CHUNK
# public rows (downloads) go to both partitions; a sample that failed before goes to stajichlab (not preemptible)
ids=()
take() {       # take NAME PART ACCT N THROTTLE < rows
  local name=$1 part=$2 acct=$3 n=$4 thr=$5 f="$STORE/lists/r${round}_$1.txt"
  head -n "$n" > "$f"; local k; k=$(wc -l < "$f"); [[ $k -gt 0 ]] || { rm -f "$f"; return 0; }
  awk -v r=$round '{print $1 "\t" r}' "$f" >> "$STORE/lists/attempts.tsv"
  ids+=("$(sbatch --parsable -p "$part" -A "$acct" $( [[ $part == preempt ]] && echo --nice=500 ) -J "sk-$SPECIES-$name" \
    --array=1-$k%$thr -c 4 --mem 8G -t 3:00:00 -o "$STORE/logs/sk.%A_%a.log" --wrap "source $STORE/env.sh; LIST=$f bash $STORE/bin/sketch_one.sh $STORE")")
  echo "round $round: $name $k samples on $part (array ${ids[-1]})"
}
awk 'FILENAME==ARGV[1]{n[$1]++; next} !($1 in n)' "$STORE/lists/attempts.tsv" "$STORE/lists/todo.now" > "$STORE/lists/fresh.now"
awk 'FILENAME==ARGV[1]{n[$1]++; next} ($1 in n)' "$STORE/lists/attempts.tsv" "$STORE/lists/todo.now" > "$STORE/lists/retry.now"
nr=$(wc -l < "$STORE/lists/retry.now"); [[ $nr -gt $rl ]] && nr=$rl
take retry stajichlab stajichlab "$nr" "${LAB_THROTTLE:-25}" < "$STORE/lists/retry.now"
rl=$(( rl - nr ))
take pre preempt preempt "$rp" "${PRE_THROTTLE:-35}" < "$STORE/lists/fresh.now"
take lab stajichlab stajichlab "$rl" "${LAB_THROTTLE:-25}" < <(tail -n +$(( rp + 1 )) "$STORE/lists/fresh.now")
if [[ ${#ids[@]} -eq 0 ]]; then
  # no room this time: look again in 30 minutes
  p=$(sbatch --parsable --begin=now+30minutes -p stajichlab -A stajichlab -c 1 --mem 2G -t 0:30:00 -J "sk-$SPECIES-plan" \
      -o "$STORE/logs/plan.%j.log" --wrap "PLAN=1 bash $STORE/bin/submit.sh $STORE")
  echo "round $round: no room under the account limits ($ntodo to do); planner again at +30 min ($p)"
  echo -e "$(date +%F_%T)\twait\t$p\t$ntodo to do" >> "$STORE/jobs.tsv"; exit 0
fi
dep=$(IFS=:; echo "afterany:${ids[*]}")
p=$(plan_job "$dep")
echo -e "$(date +%F_%T)\tround$round\t${ids[*]}\t$p\t$ntodo to do" >> "$STORE/jobs.tsv"
echo "round $round queued; next planner $p"
