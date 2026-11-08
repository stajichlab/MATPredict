#!/bin/bash
# One sample of a sketch store: a COVERAGE-x subsample of its reads, MAT-typed (MATPredict reads-type) and sketched
# (sourmash_branchwater singlesketch, k=21,k=31, scaled=100, abundances; k-mers seen once are dropped).
#
#   sketch_one.sh STORE INDEX          INDEX = data row of STORE/manifest.tsv (1-based)
#   LIST=rows.txt sketch_one.sh STORE  as an array task: the row number is line SLURM_ARRAY_TASK_ID of LIST (see submit.sh)
#
# Writes STORE/sigs/<sample_id>.sig.zip, STORE/stats/<sample_id>.tsv and, for typing, STORE/mat/<run or sample_id>.*
# Public rows: the MATPredict reads worker (STORE/bin/genotype_mat_one_accession.sh) downloads the run from SRA, takes target_spots
# spots in evenly spaced windows, screens them with GmG (SCREEN=1) and types them; its subsample is kept (FASTQ_OUT) and sketched.
# Local rows: the FASTQ files (local_pattern, shell patterns relative to LOCAL_INPUT) are read directly, the first N read pairs.
# Idempotent: a sample with a sketch and a typing result is skipped; nothing but the subsample (in node scratch) is downloaded or kept.
# Environment: PANEL_DIR, PANEL_PREFIX (the species' MATPredict reads panel), LOCAL_INPUT, GENOME_KB, COVERAGE, BW (sourmash_branchwater bin).
set -uo pipefail
STORE=${1:?STORE}; IDX=${2:-}
[[ -z $IDX && -n ${LIST:-} ]] && IDX=$(sed -n "${SLURM_ARRAY_TASK_ID:?}p" "$LIST")
[[ -n $IDX ]] || { echo "need INDEX or LIST + SLURM_ARRAY_TASK_ID" >&2; exit 2; }
BW=${BW:-/opt/linux/rocky/8.x/x86_64/pkgs/sourmash_branchwater/0.9.13/bin}
P=/bigdata/stajichlab/jstajich/projects/MATPredict
PY=${PY:-$P/.pixi/envs/default/bin/python}
export PYTHONPATH=${PYTHONPATH:-$P/.claude/worktrees/run-b317ece/src} MATPREDICT_DB_ROOT=${MATPREDICT_DB_ROOT:-$P/.claude/worktrees/run-b317ece/db}
PANEL_DIR=${PANEL_DIR:?PANEL_DIR}; PANEL_PREFIX=${PANEL_PREFIX:?PANEL_PREFIX}
GENOME_KB=${GENOME_KB:?GENOME_KB}; COVERAGE=${COVERAGE:-30}; PARAMS=${PARAMS:-k=21,k=31,scaled=100,abund}
mkdir -p "$STORE"/{sigs,stats,mat}
# fields joined with \037, not a tab: `read` merges runs of a whitespace separator, which would drop the empty fields of local rows
row=$(awk -F'\t' -v n=$((IDX + 1)) 'BEGIN{OFS="\037"} NR==1{for(i=1;i<=NF;i++)c[$i]=i} NR==n{print $c["sample_id"], $c["source"], $c["run_accession"], $c["read_count"], $c["target_spots"], $c["layout"], $c["local_pattern"], (("ena_fastq" in c) ? $c["ena_fastq"] : "")}' "$STORE/manifest.tsv")
[[ -n $row ]] || { echo "no manifest row $IDX" >&2; exit 2; }
IFS=$'\037' read -r SID SRC RUN SPOTS TSPOTS LAYOUT PAT ENA <<< "$row"
# take_pairs: the first N read pairs (N = COVERAGE x genome size) of r1[]/r2[] (gzipped FASTQ, read with $reader: cat or curl) into $fq,
# interleaved when paired
take_pairs() {
  local L N
  L=$($reader "${r1[0]}" | zcat -f 2>/dev/null | head -2 | tail -1 | tr -d '\n' | wc -c)
  [[ $L -gt 0 ]] || return 0
  if [[ ${#r2[@]} -gt 0 ]]; then
    N=$(( GENOME_KB * 1000 * COVERAGE / (2 * L) ))
    paste <(for f in "${r1[@]}"; do $reader "$f"; done | zcat -f 2>/dev/null | head -n $((4 * N)) | paste - - - -) \
          <(for f in "${r2[@]}"; do $reader "$f"; done | zcat -f 2>/dev/null | head -n $((4 * N)) | paste - - - -) | tr '\t' '\n' > "$fq"
  else
    N=$(( GENOME_KB * 1000 * COVERAGE / L ))
    for f in "${r1[@]}"; do $reader "$f"; done | zcat -f 2>/dev/null | head -n $((4 * N)) > "$fq"
  fi
  # a stream that broke mid-file (seen on ENA HTTPS for 3 GB 2013 files) leaves a truncated record that sourmash rejects:
  # keep only whole records (4 lines, '@' header, sequence and quality of equal length), in whole pairs when interleaved
  awk -v per=$(( ${#r2[@]} > 0 ? 2 : 1 )) 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0}
    NR%4==0{ if (substr(h,1,1)!="@" || length(s)!=length($0)) exit; b=b h "\n" s "\n" p "\n" $0 "\n"; if (++k==per) {printf "%s", b; b=""; k=0} }' \
    "$fq" > "$fq.ok" && mv "$fq.ok" "$fq"
}
KEY=${RUN:-$SID}
typed() { [[ -s $STORE/mat/$KEY.$PANEL_PREFIX.tsv ]]; }
[[ -s $STORE/sigs/$SID.sig.zip ]] && typed && { echo "skip $SID"; exit 0; }
work=$(mktemp -d "${SCRATCH:-/tmp}/sk.$SID.XXXXXX"); trap 'rm -rf "$work"' EXIT
fq=$work/sub.fq; t0=$(date +%s); state=ok

if [[ $SRC == public ]]; then
  # the worker types into STORE/mat; if the run is typed already (an earlier attempt), it only extracts (SKIP_TYPE=1, scratch OUT)
  out=$STORE/mat; extra=(); typed && { out=$work/retype; extra=(SKIP_TYPE=1); }
  env "${extra[@]}" TARGET_SPOTS="$TSPOTS" SCREEN=${SCREEN:-1} SKETCH=0 FASTQ_OUT="$fq" FALLBACK_K=${FALLBACK_K:-21} \
    bash "$STORE/bin/genotype_mat_one_accession.sh" "$RUN" "$SPOTS" "$PANEL_DIR" "$PANEL_PREFIX" "$out" aspergillus_fumigatus "$LAYOUT" \
    > "$work/worker.log" 2>&1 || state=worker_failed
  if [[ ! -s $fq && -n $ENA ]]; then
    # not at NCBI/AWS (recent ENA submissions are mirrored late): stream the run's FASTQ from ENA over HTTPS and keep the first N pairs
    # (only those bytes are downloaded). No GmG screen on this path; the sketch QC (query qc) still flags mixtures. Typed below.
    IFS=';' read -ra urls <<< "$ENA"; r1=("https://${urls[0]}"); r2=(); [[ ${#urls[@]} -ge 2 ]] && r2=("https://${urls[-1]}")
    rm -f "$fq"; reader="curl -sfL --retry 3 --retry-delay 20"; take_pairs; ENA_USED=1
    # a broken stream that delivered under half the target (and under half the run, for runs smaller than the target) is retried by
    # the planner rather than sketched thin
    got=$(awk 'NR%4==2{s+=length($0)}END{print s+0}' "$fq" 2>/dev/null || echo 0)
    want=$(awk -F'\t' -v n=$((IDX + 1)) -v t=$(( GENOME_KB * 1000 * COVERAGE )) 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} NR==n{b=$c["base_count"]+0; print (b && b < t) ? b : t}' "$STORE/manifest.tsv")
    if [[ -s $fq && $got -lt $(( want / 2 )) ]]; then
      echo "$SID: ENA stream gave $got of $want bases; leaving it for a retry" >&2
      echo -e "$SID\t$SRC\t$KEY\t$got\t0\t0\tena_partial" > "$STORE/stats/$SID.tsv"; exit 1
    fi
    if [[ -s $fq ]]; then state=ok; echo -e "{\"run\":\"$RUN\",\"state\":\"done_ena\",\"source\":\"ena-https\"}" > "$STORE/mat/$KEY.ena.json"; fi
  fi
  [[ -s $fq ]] || { echo "$SID: no subsample ($state): $(tail -3 "$work/worker.log" | tr '\n' ' ')" >&2; \
    echo -e "$SID\t$SRC\t$KEY\tno_reads\t0\t0\t$state" > "$STORE/stats/$SID.tsv"; exit 1; }
else
  IN=${LOCAL_INPUT:?LOCAL_INPUT}; r1=(); r2=()
  IFS=';' read -ra parts <<< "$PAT"
  for p in "${parts[@]}"; do ex=( $IN/$p ); r1+=("${ex[0]}"); [[ ${#ex[@]} -ge 2 ]] && r2+=("${ex[1]}"); done
  for f in "${r1[@]}" "${r2[@]}"; do [[ -s $f ]] || { echo "$SID: missing $f" >&2; echo -e "$SID\t$SRC\t$KEY\tmissing_files\t0\t0\tmissing" > "$STORE/stats/$SID.tsv"; exit 1; }; done
  reader="cat"; take_pairs
fi
if [[ $SRC != public || ${ENA_USED:-0} == 1 ]]; then
  if ! typed; then
    printf 'sample\treads\n%s\t%s\n' "$SID" "$fq" > "$work/samples.tsv"
    idio=(); for n in ${IDIOMORPHS:-MAT1-1 MAT1-2}; do idio+=(--idiomorph "$n=$PANEL_DIR/${LOCAL_LINEAGE:-L1}/${PANEL_PREFIX}_$n.fasta"); done
    $PY -m MATPredict reads-type "${idio[@]}" --samples "$work/samples.tsv" --out "$STORE/mat/$KEY.$PANEL_PREFIX.tsv" > "$work/type.log" 2>&1 || { state=type_failed; rm -f "$STORE/mat/$KEY.$PANEL_PREFIX.tsv"; }
  fi
fi
t1=$(date +%s)
bases=$(awk 'NR%4==2{s+=length($0)}END{print s+0}' "$fq")
if [[ ! -s $STORE/sigs/$SID.sig.zip ]]; then
  "$BW/sourmash" scripts singlesketch -q -p "$PARAMS" -n "$SID" -o "$work/raw.sig.zip" "$fq" > "$work/sketch.log" 2>&1 \
    && "$BW/sourmash" sig filter -q -m 2 "$work/raw.sig.zip" -o "$work/f.sig.zip" >> "$work/sketch.log" 2>&1 \
    && mv "$work/f.sig.zip" "$STORE/sigs/$SID.sig.zip" || { state=sketch_failed; echo "$SID sketch failed: $(tail -2 "$work/sketch.log")" >&2; }
fi
raw=$(stat -c %s "$work/raw.sig.zip" 2>/dev/null || echo 0)
echo -e "$SID\t$SRC\t$KEY\t$bases\t$(awk -v b=$bases -v g=$GENOME_KB 'BEGIN{printf "%.1f", b/(g*1000)}')\t$raw\t$state\t$((t1 - t0))\t$(( $(date +%s) - t1 ))" > "$STORE/stats/$SID.tsv"
echo "$SID $state: $((t1 - t0)) s reads, $(( $(date +%s) - t1 )) s sketch, $bases bases"
[[ $state == ok ]]
