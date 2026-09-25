#!/usr/bin/bash -l
#SBATCH --job-name=matpredict-panel
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --output=logs/panel_%A.log
#
# Run `matpredict detect` over every BFD genome in one clade.
#
# The clade is chosen by any column of samples.csv, so the same script covers a
# genus, a family, an order or a class:
#
#   cd /bigdata/stajichlab/jstajich/projects/MATPredict
#   sbatch scripts/run_clade_panel.slurm                      # Tremellales (default)
#   sbatch --export=ALL,RANK=GENUS,CLADE=Saccharomyces        scripts/run_clade_panel.slurm
#   sbatch --export=ALL,RANK=GENUS,CLADE=Zygosaccharomyces    scripts/run_clade_panel.slurm
#   sbatch --export=ALL,RANK=CLASS,CLADE=Saccharomycetes,PARTITION=stajichlab \
#          --partition=stajichlab --time=12:00:00             scripts/run_clade_panel.slurm
#
# NOTE on "the whole Saccharomycotina subphylum": samples.csv has a SUBPHYLUM
# column but it is EMPTY in all 23,683 rows. Saccharomycotina is split across
# CLASS values there -- Saccharomycetes (2,648) and Pichiomycetes (2,967) --
# so select by CLASS, or run it twice. Note also that `MATsc` is scoped to
# taxid 4893 = Saccharomycetaceae, so genomes outside that family fall to
# phylum_fallback and are searched against every Ascomycota family.
#
# Submit from the repo root: SLURM opens logs/tremellales_%A.log relative to
# the submission directory and does so BEFORE this script runs, so the job
# fails immediately if logs/ is missing. It is committed with a .gitkeep.
#
# Override anything on the command line, e.g. a different clade or output dir:

# Measured budget: ~20 s/genome on a chromosome-level assembly after the
# 2026-09-21 roster trim (STE20 off the search list took a run from 535 s to
# 19 s). 343 Tremellales genomes at 8-way parallelism is roughly 15-45 min,
# which fits `short` (2 h) with headroom. If you widen ORDER to a clade with
# thousands of genomes, move to --partition=stajichlab and raise --time.
#
# NOTE: this script never uses ${BASH_SOURCE[0]} -- it resolves to the wrong
# path inside a SLURM work directory. REPO is an absolute path, overridable.
#
# ALSO: do not edit the working tree while this is running. It executes
# `$REPO/src` directly, so a half-saved file becomes a SyntaxError in every
# concurrently-running genome. That cost 277 of 2,283 genomes on 2026-09-21.
# Use a git worktree for edits, or wait for the job.

set -uo pipefail

REPO="${REPO:-/bigdata/stajichlab/jstajich/projects/MATPredict}"
SAMPLES="${SAMPLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv}"
LIBRARY="${LIBRARY:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes}"
RANK="${RANK:-ORDER}"          # any column of samples.csv: GENUS, FAMILY, ORDER, CLASS...
CLADE="${CLADE:-${ORDER:-Tremellales}}"
OUT_DIR="${OUT_DIR:-$REPO/results/$(date +%Y-%m-%d)_${CLADE,,}_panel}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

PY="$REPO/.pixi/envs/default/bin/python"
export PATH="$REPO/.pixi/envs/default/bin:$PATH"   # exonerate, miniprot, tblastn
# SRC lets a run execute a WORKTREE's source while still using $REPO's pixi
# environment and results directory -- how a change is panel-tested without
# touching the live tree (the 2026-09-21 lesson).
SRC="${SRC:-$REPO/src}"
# FAIL LOUDLY when SRC is not visible from the compute node. A $PYTHONPATH
# pointing at a missing directory does NOT error: python silently falls back to
# the INSTALLED MATPredict, so the job runs to completion against the WRONG
# code and its output looks entirely normal. That cost two launches on
# 2026-09-23, when a worktree on /scratch -- which is NODE-LOCAL -- existed
# only on the submitting node. A worktree for a SLURM run belongs on /bigdata.
if [[ ! -f "$SRC/MATPredict/__init__.py" ]]; then
    echo "SRC=$SRC is not a MATPredict source tree on this node ($(hostname))" >&2
    echo "  /scratch is node-local; put a worktree for a SLURM run on /bigdata" >&2
    exit 1
fi
export PYTHONPATH="$SRC"
# The curated database must come from the SAME tree as the code. `detect`
# reads $MATPREDICT_DB_ROOT, else <working dir>/db -- and a SLURM job's working
# dir is wherever it was submitted from, usually the main checkout. That made a
# 2026-09-25 re-run search the main checkout's OLD database while running the
# worktree's code, so two newly curated records never reached the search and
# the run reported no change. Pinned next to SRC unless set explicitly.
export MATPREDICT_DB_ROOT="${MATPREDICT_DB_ROOT:-$(dirname "$SRC")/db}"
if [[ ! -f "$MATPREDICT_DB_ROOT/Ascomycota/order.yml" ]]; then
    echo "MATPREDICT_DB_ROOT=$MATPREDICT_DB_ROOT is not a curated database on $(hostname)" >&2
    exit 1
fi

# $SCRATCH is node-local NVMe, pre-created per job, and NOT shared between
# nodes. Genomes are decompressed here and deleted straight after use; only the
# small YAML reports are copied back to shared storage.
WORK="${SCRATCH:?SCRATCH is not set -- run this under SLURM}/matpredict_${CLADE//[^A-Za-z0-9]/_}_${SLURM_JOB_ID:-local}"

mkdir -p "$WORK/fasta" "$OUT_DIR/runs" "$REPO/logs"

echo "repo=$REPO"
echo "src=$SRC"
echo "db=$MATPREDICT_DB_ROOT"
echo "clade=$RANK=$CLADE  samples=$SAMPLES"
echo "scratch=$WORK"
echo "out=$OUT_DIR  threads=$THREADS"

# ---------------------------------------------------------------- sample list
# One "ASMID<TAB>TAXID" line per genome that exists on disk.
#
# SAMPLE_LIST, when set, is used INSTEAD of the RANK/CLADE selection: a
# prepared "ASMID<TAB>TAXID" file (on /bigdata), e.g. a stratified pilot pick.
# Name the run with CLADE anyway, since it names the scratch dir and OUT_DIR.
if [[ -n "${SAMPLE_LIST:-}" ]]; then
    [[ -s "$SAMPLE_LIST" ]] || { echo "SAMPLE_LIST=$SAMPLE_LIST missing or empty on $(hostname)" >&2; exit 1; }
    while IFS=$'\t' read -r asmid taxid; do
        [[ -n "$taxid" && -e "$LIBRARY/$asmid.fa.gz" ]] && printf '%s\t%s\n' "$asmid" "$taxid"
    done < "$SAMPLE_LIST" > "$WORK/panel.tsv"
    echo "SAMPLE_LIST=$SAMPLE_LIST: $(wc -l < "$SAMPLE_LIST") rows, $(wc -l < "$WORK/panel.tsv") with a genome file and taxid" >&2
else
"$PY" - "$SAMPLES" "$LIBRARY" "$RANK" "$CLADE" > "$WORK/panel.tsv" <<'PYEOF'
import csv, os, sys
samples, library, rank, clade = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
n = kept = 0
reader = csv.DictReader(open(samples))
if rank not in (reader.fieldnames or []):
    sys.exit(f"no column {rank!r} in {samples}; have: {reader.fieldnames}")
for r in reader:
    if (r.get(rank) or "") != clade:
        continue
    n += 1
    asmid, taxid = r["ASMID"], (r.get("NCBI_TAXONID") or "").strip()
    if not taxid:
        continue                       # no taxid: detect could not route it
    if os.path.exists(os.path.join(library, asmid + ".fa.gz")):
        print(f"{asmid}\t{taxid}")
        kept += 1
print(f"{rank}={clade}: {n} rows, {kept} with a genome file", file=sys.stderr)
PYEOF
fi

# ----------------------------------------------------------- taxonomy warm-up
# Resolve every distinct taxid ONCE, serially, before fanning out. The routing,
# phylum-fallback and genetic-code lookups all read one cached efetch document
# per taxid, so afterwards the parallel runs are cache hits. Without this, 128
# concurrent processes shared NCBI's keyless 3 req/s; failed lookups silently
# widened routing (one taxid routed three different ways in one pilot, several
# genomes ran `exhaustive` at ~30 min each). The cache dir is pinned so every
# job shares it regardless of where it was submitted from.
export MATPREDICT_CACHE_DIR="${MATPREDICT_CACHE_DIR:-$REPO/.matpredict_cache}"
cut -f2 "$WORK/panel.tsv" | sort -u | "$PY" -c '
import sys, time
from MATPredict.db.taxonomy import (default_genetic_code, default_lineage_phylum_name,
                                    default_lineage_taxids)
bad = 0
for line in sys.stdin:
    t = int(line)
    try:
        default_lineage_taxids(t); default_lineage_phylum_name(t); default_genetic_code(t)
    except Exception as exc:
        bad += 1
        print(f"taxonomy warm-up: taxid {t} failed: {exc}", file=sys.stderr)
    time.sleep(0.35)
print(f"taxonomy warm-up done, {bad} failed", file=sys.stderr)
'

TOTAL=$(wc -l < "$WORK/panel.tsv")
if [[ "$TOTAL" -eq 0 ]]; then
    echo "no genomes matched $RANK=$CLADE -- check the spelling against samples.csv" >&2
    exit 1
fi
echo "genomes to run: $TOTAL"

# ------------------------------------------------------------------ one genome
run_one() {
    asmid="$1"; taxid="$2"
    out="$OUT_DIR/runs/$asmid"
    # Resumable: a finished genome is skipped, so a timed-out job can simply be
    # resubmitted with the same OUT_DIR.
    [[ -s "$out/detection_report.yaml" ]] && return 0

    fna="$WORK/fasta/$asmid.fna"
    if ! zcat "$LIBRARY/$asmid.fa.gz" > "$fna" 2>"$WORK/$asmid.unzip.err"; then
        echo "FAIL $asmid decompress" >&2; rm -f "$fna"; return 0
    fi

    mkdir -p "$out"
    # NOT `if ! ...; then echo rc=$?`: inside that block `$?` is the status of
    # the negated test, which is always 0, so every failure reported rc=0 and
    # the real exit code was lost. Capture it from the command itself.
    # Wall time per genome, written beside the report. Without it runtime had
    # to be reconstructed from file mtimes, which a non-empty stderr breaks.
    t0=$(date +%s)
    timeout 3600 "$PY" -m MATPredict detect \
        --genome "$fna" --taxid "$taxid" --out-dir "$out" \
        --evidence-diagnostics "$out/evidence_diagnostics.jsonl" \
        > "$out/stdout.log" 2> "$out/stderr.log"
    rc=$?
    echo "$(( $(date +%s) - t0 ))" > "$out/wall_seconds"
    if [[ $rc -ne 0 ]]; then
        # The last stderr line is usually the whole story (a SyntaxError, a
        # missing tool, an OOM); printing it here saves opening 277 files.
        echo "FAIL $asmid detect rc=$rc: $(tail -1 "$out/stderr.log" 2>/dev/null)" >&2
    fi
    # Node-local scratch is not unlimited and these are ~20-200 MB each.
    rm -f "$fna"
}
export -f run_one
export REPO LIBRARY OUT_DIR WORK PY

# xargs, not a job array: the per-genome cost is ~20 s, far below the point
# where a separate SLURM job per unit would pay for its own scheduling.
awk -F'\t' '{print $1" "$2}' "$WORK/panel.tsv" \
    | xargs -P "$THREADS" -L1 bash -c 'run_one "$0" "$1"'

# --------------------------------------------------------------------- summary
DONE=$(find "$OUT_DIR/runs" -name detection_report.yaml -size +0 | wc -l)
echo "completed $DONE / $TOTAL"

"$PY" -m MATPredict detect rollout-summary \
    --reports-dir "$OUT_DIR/runs" --out "$OUT_DIR/rollout_summary.yaml" \
    2>/dev/null || echo "(rollout-summary skipped; run it by hand if you want the tallies)"

# Reports are small; keep them compressed on shared storage.
tar cf - -C "$OUT_DIR/runs" . 2>/dev/null \
    | zstd -T0 -q -o "$OUT_DIR/reports.tar.zst" \
    && echo "archived $OUT_DIR/reports.tar.zst"

rm -rf "$WORK/fasta"
echo "done."
