#!/usr/bin/env bash
# Genotype the MAT locus of ONE SRA run (one accession): which idiomorph does it carry? Matches its reads against a species panel. The unit of work for every runner: the Slurm wrapper (stream_type_sra.slurm), submit_species.sh, or a
# Nextflow process. It knows nothing about Slurm: everything comes in as arguments or environment, everything goes out as two small files.
#
# Usage: genotype_mat_one_accession.sh RUN SPOTS PANEL_DIR PANEL_PREFIX OUT_DIR [GROUP] [LAYOUT]
#   RUN           SRA/ENA/DDBJ run accession (SRR, ERR, DRR)
#   SPOTS         number of spots (ENA/SRA read_count); 0 or empty = read the whole run
#   PANEL_DIR     folder with PANEL_PREFIX_<IDIOMORPH>.fasta (IDIOMORPHS="MAT1-1 MAT1-2" by default)
#   OUT_DIR       receives <run>.<prefix>.tsv (the reads-type row) and <run>.run.json (the run record)
# Environment (all optional): TARGET_SPOTS=2000000  WINDOWS=8  SRA_URL=https://sra-pub-run-odp.s3.amazonaws.com/sra  IDIOMORPHS
#   SCRATCH=<fast local dir for the temporary .sra and FASTQ, default /tmp>  KEEP_SRA=1 keeps the temporary files
#   FASTQ_OUT=path keeps the subsample FASTQ there; SKIP_TYPE=1 stops before typing (follow-up analyses). MAX_SRA_SIZE=50G caps the prefetch fallback. SRATOOLKIT_BIN, PY, PYTHONPATH, MATPREDICT_DB_ROOT: where fastq-dump and `python -m MATPredict` come from (defaults: the HPCC installs).
# Exit: 0 done (or already done: the result file is the done-marker), 1 failed (the state is in the run.json), 2 bad arguments. Idempotent; safe to retry.
# Not kept: the .sra and every FASTQ. Cost per run measured 2026-10-26: about 1 GB download, 110 s, 5 GB RAM, up to 4 CPUs (the typing step).
set -uo pipefail
RUN=${1:?RUN}; SPOTS=${2:-0}; PDIR=${3:?PANEL_DIR}; PFX=${4:?PANEL_PREFIX}; OUT=${5:?OUT_DIR}; GROUP=${6:-}; LAYOUT=${7:-}
TARGET=${TARGET_SPOTS:-2000000}; WIN=${WINDOWS:-8}; SRA_URL=${SRA_URL:-https://sra-pub-run-odp.s3.amazonaws.com/sra}
P=/bigdata/stajichlab/jstajich/projects/MATPredict
export PATH=${SRATOOLKIT_BIN:-/opt/linux/rocky/8.x/x86_64/pkgs/sratoolkit/3.2.0/bin}:$PATH
PY=${PY:-$P/.pixi/envs/default/bin/python}
export PYTHONPATH=${PYTHONPATH:-$P/.claude/worktrees/run-b317ece/src}
export MATPREDICT_DB_ROOT=${MATPREDICT_DB_ROOT:-$P/.claude/worktrees/run-b317ece/db}
mkdir -p "$OUT"
[[ -s $OUT/$RUN.$PFX.tsv ]] && { echo "skip $RUN"; exit 0; }
work=$(mktemp -d "${SCRATCH:-/tmp}/sra_type.$RUN.XXXXXX")
state=started; t0=$(date +%s)
json() { printf '{"run":"%s","group":"%s","layout":"%s","spots_total":%s,"target_spots":%s,"windows":%s,"sra_bytes":%s,"secs_download":%s,"secs_convert":%s,"secs_type":%s,"state":"%s","source":"%s","host":"%s"}\n' \
  "$RUN" "$GROUP" "$LAYOUT" "${SPOTS:-0}" "$TARGET" "$WIN" "${bytes:-0}" "${d1:-0}" "${c1:-0}" "${y1:-0}" "$state" "${src:-none}" "$(hostname)" > "$OUT/$RUN.run.json"; }
trap 'json; [[ ${KEEP_SRA:-0} == 1 ]] || rm -rf "$work"' EXIT
# 1. download (resumable, 3 tries) from NCBI's copy on AWS; a 404 means that run is not mirrored there (some European runs are not): then ask NCBI itself through the SRA toolkit
ok=0; src=aws-open-data
for a in 1 2 3; do curl -sS --fail -L -C - -o "$work/$RUN.sra" "$SRA_URL/$RUN/$RUN" && { ok=1; break; }; sleep $((a*10)); done
if [[ $ok != 1 ]]; then
  rm -f "$work/$RUN.sra"; src=ncbi-prefetch
  if prefetch --max-size "${MAX_SRA_SIZE:-50G}" -O "$work" "$RUN" > "$work/prefetch.log" 2>&1 && [[ -s $work/$RUN/$RUN.sra ]]; then mv "$work/$RUN/$RUN.sra" "$work/$RUN.sra"; ok=1; fi
fi
# neither AWS nor NCBI has it (an ENA-only run, e.g. ERR2188722 on 2026-10-26: NCBI answers "no data (404)"): say so, do not call it a failure of the pipeline
if [[ $ok != 1 ]] && grep -q "no data" "$work/prefetch.log" 2>/dev/null; then state=unavailable_at_ncbi; fi
d1=$(( $(date +%s) - t0 )); bytes=$(stat -c %s "$work/$RUN.sra" 2>/dev/null || echo 0)
[[ $ok == 1 && $bytes -gt 0 ]] || { [[ $state == unavailable_at_ncbi ]] || state=failed_download; exit 1; }
# 2. evenly spaced windows of spots (1-based, inclusive); a run no larger than the target is read whole
t1=$(date +%s); fq="$work/sub.fastq"; : > "$fq"
if [[ ${SPOTS:-0} -le $TARGET ]]; then
  fastq-dump -Z --split-spot --ncbi_error_report never "$work/$RUN.sra" 2>/dev/null >> "$fq"
else
  w=$(( TARGET / WIN ))
  for i in $(seq 0 $((WIN-1))); do
    s=$(( (2*i+1) * SPOTS / (2*WIN) - w/2 + 1 )); [[ $s -lt 1 ]] && s=1; e=$(( s + w - 1 ))
    fastq-dump -N $s -X $e -Z --split-spot --ncbi_error_report never "$work/$RUN.sra" 2>/dev/null >> "$fq"
  done
fi
c1=$(( $(date +%s) - t1 )); nlines=$(wc -l < "$fq")
[[ $nlines -gt 0 ]] || { state=failed_convert; exit 1; }
rm -f "$work/$RUN.sra"
# optional, for follow-up analyses: keep a copy of the subsample (FASTQ_OUT=path) and/or stop before typing (SKIP_TYPE=1)
[[ -n ${FASTQ_OUT:-} ]] && cp "$fq" "$FASTQ_OUT"
[[ ${SKIP_TYPE:-0} == 1 ]] && { state=done_untyped; exit 0; }
# 3. type
t2=$(date +%s)
printf 'sample\treads\n%s\t%s\n' "$RUN" "$fq" > "$work/samples.tsv"
idio=(); for n in ${IDIOMORPHS:-MAT1-1 MAT1-2}; do idio+=(--idiomorph "$n=$PDIR/${PFX}_$n.fasta"); done
$PY -m MATPredict reads-type "${idio[@]}" --samples "$work/samples.tsv" --out "$work/$RUN.$PFX.tsv" || { state=failed_type; y1=$(( $(date +%s) - t2 )); exit 1; }
y1=$(( $(date +%s) - t2 )); state=done
mv "$work/$RUN.$PFX.tsv" "$OUT/$RUN.$PFX.tsv"
echo "done $RUN: $((nlines/4)) reads typed, download ${d1}s ($((bytes/1000000)) MB), convert ${c1}s, type ${y1}s"
