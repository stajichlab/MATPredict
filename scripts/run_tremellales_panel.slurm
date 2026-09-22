#!/usr/bin/bash -l
#SBATCH --job-name=matpredict-tremellales
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=2:00:00
#SBATCH --output=logs/tremellales_%A.log
#
# Run `matpredict detect` over every Tremellales genome in the BFD library.
#
#   cd /bigdata/stajichlab/jstajich/projects/MATPredict
#   sbatch scripts/run_tremellales_panel.slurm
#
# Submit from the repo root: SLURM opens logs/tremellales_%A.log relative to
# the submission directory and does so BEFORE this script runs, so the job
# fails immediately if logs/ is missing. It is committed with a .gitkeep.
#
# Override anything on the command line, e.g. a different clade or output dir:
#   sbatch --export=ALL,ORDER=Agaricales,OUT_DIR=$PWD/results/agaricales \
#          scripts/run_tremellales_panel.slurm
#
# Measured budget: ~20 s/genome on a chromosome-level assembly after the
# 2026-09-21 roster trim (STE20 off the search list took a run from 535 s to
# 19 s). 343 Tremellales genomes at 8-way parallelism is roughly 15-45 min,
# which fits `short` (2 h) with headroom. If you widen ORDER to a clade with
# thousands of genomes, move to --partition=stajichlab and raise --time.
#
# NOTE: this script never uses ${BASH_SOURCE[0]} -- it resolves to the wrong
# path inside a SLURM work directory. REPO is an absolute path, overridable.

set -uo pipefail

REPO="${REPO:-/bigdata/stajichlab/jstajich/projects/MATPredict}"
SAMPLES="${SAMPLES:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv}"
LIBRARY="${LIBRARY:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes}"
ORDER="${ORDER:-Tremellales}"
OUT_DIR="${OUT_DIR:-$REPO/results/$(date +%Y-%m-%d)_${ORDER,,}_panel}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

PY="$REPO/.pixi/envs/default/bin/python"
export PATH="$REPO/.pixi/envs/default/bin:$PATH"   # exonerate, miniprot, tblastn
export PYTHONPATH="$REPO/src"

# $SCRATCH is node-local NVMe, pre-created per job, and NOT shared between
# nodes. Genomes are decompressed here and deleted straight after use; only the
# small YAML reports are copied back to shared storage.
WORK="${SCRATCH:?SCRATCH is not set -- run this under SLURM}/matpredict_${ORDER}_${SLURM_JOB_ID:-local}"

mkdir -p "$WORK/fasta" "$OUT_DIR/runs" "$REPO/logs"

echo "repo=$REPO"
echo "order=$ORDER  samples=$SAMPLES"
echo "scratch=$WORK"
echo "out=$OUT_DIR  threads=$THREADS"

# ---------------------------------------------------------------- sample list
# One "ASMID<TAB>TAXID" line per genome that exists on disk.
"$PY" - "$SAMPLES" "$LIBRARY" "$ORDER" > "$WORK/panel.tsv" <<'PYEOF'
import csv, os, sys
samples, library, order = sys.argv[1], sys.argv[2], sys.argv[3]
n = kept = 0
for r in csv.DictReader(open(samples)):
    if (r.get("ORDER") or "") != order:
        continue
    n += 1
    asmid, taxid = r["ASMID"], (r.get("NCBI_TAXONID") or "").strip()
    if not taxid:
        continue                       # no taxid: detect could not route it
    if os.path.exists(os.path.join(library, asmid + ".fa.gz")):
        print(f"{asmid}\t{taxid}")
        kept += 1
print(f"{order}: {n} rows, {kept} with a genome file", file=sys.stderr)
PYEOF

TOTAL=$(wc -l < "$WORK/panel.tsv")
if [[ "$TOTAL" -eq 0 ]]; then
    echo "no genomes matched ORDER=$ORDER -- check the spelling against samples.csv" >&2
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
    if ! timeout 3600 "$PY" -m MATPredict detect \
            --genome "$fna" --taxid "$taxid" --out-dir "$out" \
            --evidence-diagnostics "$out/evidence_diagnostics.jsonl" \
            > "$out/stdout.log" 2> "$out/stderr.log"; then
        echo "FAIL $asmid detect rc=$?" >&2
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
