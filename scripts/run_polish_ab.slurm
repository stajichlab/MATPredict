#!/usr/bin/bash -l
#SBATCH --job-name=matpredict-polish-ab
#SBATCH --partition=stajichlab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=6:00:00
#SBATCH --output=logs/polish_ab_%A.log
#
# A/B the polish idiomorph-scope cut (cut 1) on a stratified panel sample.
#
# Each genome is decompressed ONCE and run through BOTH arms against that same
# file, so the comparison is not confounded by assembly differences:
#
#   BASE  a detached worktree pinned to the commit the cut was made against
#   CUT1  the worktree carrying the cut
#
# Both arms are worktrees, never $REPO/src, so this job is safe to run while
# other panels execute against the main tree -- the lesson of 2026-09-21, when
# editing a live tree killed 277 of 2,283 genomes.
#
# ARM ORDER ALTERNATES per genome. Running one arm first warms the OS page
# cache for the genome FASTA, which is a real advantage of a few seconds for
# whichever arm runs second. Alternating cancels it across the sample instead
# of handing it to one arm 60 times.
#
# Timing is wall-clock on a SHARED node, so absolute seconds are not portable.
# The paired difference is, because both arms of a pair run back-to-back under
# the same load.
#
# Submit from the repo root (SLURM opens logs/ before the script runs):
#   sbatch --export=ALL,SAMPLE=/path/to/ab_sample.tsv scripts/run_polish_ab.slurm
#
# NOTE: never ${BASH_SOURCE[0]} -- it resolves wrongly in a SLURM work dir.

set -uo pipefail

REPO="${REPO:-/bigdata/stajichlab/jstajich/projects/MATPredict}"
BASE_SRC="${BASE_SRC:-/scratch/jstajich/28984668/wt-baseline/src}"
CUT1_SRC="${CUT1_SRC:-/scratch/jstajich/28984668/wt-polish/src}"
SAMPLE="${SAMPLE:-/scratch/jstajich/28984668/ab_sample.tsv}"
LIBRARY="${LIBRARY:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes}"
OUT_DIR="${OUT_DIR:-$REPO/results/$(date +%Y-%m-%d)_polish_ab}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

PY="$REPO/.pixi/envs/default/bin/python"
export PATH="$REPO/.pixi/envs/default/bin:$PATH"   # exonerate, miniprot, tblastn, makeblastdb

WORK="${SCRATCH:?SCRATCH is not set -- run this under SLURM}/matpredict_ab_${SLURM_JOB_ID:-local}"
mkdir -p "$WORK/fasta" "$OUT_DIR/runs" "$REPO/logs"

echo "repo=$REPO"
echo "base=$BASE_SRC"
echo "cut1=$CUT1_SRC"
echo "sample=$SAMPLE  out=$OUT_DIR  threads=$THREADS  scratch=$WORK"

TOTAL=$(wc -l < "$SAMPLE")
echo "genomes to run: $TOTAL"

run_pair() {
    asmid="$1"; taxid="$2"; idx="$3"
    row="$OUT_DIR/runs/$asmid/pair.tsv"
    # Resumable, like run_clade_panel.slurm: a finished pair is skipped, so a
    # timed-out job can be resubmitted with the same OUT_DIR.
    [[ -s "$row" ]] && return 0

    fna="$WORK/fasta/$asmid.fna"
    if ! zcat "$LIBRARY/$asmid.fa.gz" > "$fna" 2>"$WORK/$asmid.unzip.err"; then
        echo "FAIL $asmid decompress" >&2; rm -f "$fna"; return 0
    fi
    mkdir -p "$OUT_DIR/runs/$asmid"

    one_arm() {   # arm_name src_dir -> prints elapsed seconds, or "NA"
        local arm="$1" src="$2" out="$OUT_DIR/runs/$asmid/$1" t0 t1 rc
        mkdir -p "$out"
        t0=$(date +%s.%N)
        PYTHONPATH="$src" timeout 3600 "$PY" -m MATPredict detect \
            --genome "$fna" --taxid "$taxid" --out-dir "$out" \
            --evidence-diagnostics "$out/evidence_diagnostics.jsonl" \
            > "$out/stdout.log" 2> "$out/stderr.log"
        rc=$?
        t1=$(date +%s.%N)
        if [[ $rc -ne 0 ]]; then
            echo "FAIL $asmid $arm rc=$rc: $(tail -1 "$out/stderr.log" 2>/dev/null)" >&2
            echo "NA"; return 0
        fi
        awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}'
    }

    # Alternate which arm pays the cold-cache cost.
    if (( idx % 2 == 0 )); then
        base_s=$(one_arm base "$BASE_SRC"); cut1_s=$(one_arm cut1 "$CUT1_SRC"); order="base_first"
    else
        cut1_s=$(one_arm cut1 "$CUT1_SRC"); base_s=$(one_arm base "$BASE_SRC"); order="cut1_first"
    fi

    b="$OUT_DIR/runs/$asmid/base"; c="$OUT_DIR/runs/$asmid/cut1"
    same_report=$(cmp -s "$b/detection_report.yaml" "$c/detection_report.yaml" && echo yes || echo no)
    same_gff=$(cmp -s "$b/detected_loci.gff3" "$c/detected_loci.gff3" && echo yes || echo no)
    # NOT `grep -c ... || echo 0`: grep exits 1 when the count is 0, having
    # ALREADY printed "0", so the fallback appends a second "0" and the row
    # gains a stray line. That split 2 of 60 rows on the first run and crashed
    # the summary with an IndexError. `|| true` keeps the single printed count.
    narrowed=$(grep -c '"kind": "polish_scope"' "$c/evidence_diagnostics.jsonl" 2>/dev/null || true)
    narrowed=${narrowed:-0}
    skipped=$("$PY" - "$c/evidence_diagnostics.jsonl" <<'PYEOF'
import json, sys
n = 0
try:
    for line in open(sys.argv[1]):
        d = json.loads(line)
        if d.get("kind") == "polish_scope":
            n += len(d["rescues_skipped"])
except Exception:
    pass
print(n)
PYEOF
)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$asmid" "$taxid" "$base_s" "$cut1_s" "$same_report" "$same_gff" "$narrowed" "$skipped" > "$row"
    rm -f "$fna"
}
export -f run_pair
export REPO BASE_SRC CUT1_SRC LIBRARY OUT_DIR WORK PY

awk -F'\t' '{print $1" "$2" "NR}' "$SAMPLE" \
    | xargs -P "$THREADS" -L1 bash -c 'run_pair "$0" "$1" "$2"'

# ------------------------------------------------------------------ summary
SUM="$OUT_DIR/pairs.tsv"
{ printf 'asmid\ttaxid\tbase_s\tcut1_s\tsame_report\tsame_gff\tclusters_narrowed\tpairs_skipped\n'
  cat "$OUT_DIR"/runs/*/pair.tsv 2>/dev/null
} > "$SUM"
echo "wrote $SUM"

"$PY" - "$SUM" <<'PYEOF'
import statistics, sys
rows = [f for f in (l.rstrip("\n").split("\t") for l in open(sys.argv[1]))
        if len(f) == 8 and f[0] != "asmid"]
ok = [r for r in rows if r[2] != "NA" and r[3] != "NA"]
if not ok:
    sys.exit("no completed pairs")
b = [float(r[2]) for r in ok]; c = [float(r[3]) for r in ok]
diff = [x - y for x, y in zip(b, c)]
changed = [r[0] for r in ok if r[4] == "no" or r[5] == "no"]
print(f"pairs completed      : {len(ok)} of {len(rows)}")
print(f"baseline total       : {sum(b):.1f} s")
print(f"cut1 total           : {sum(c):.1f} s")
print(f"saved                : {sum(diff):.1f} s  ({100*sum(diff)/sum(b):.1f}%)")
print(f"per-genome saving    : median {statistics.median(diff):.1f} s, "
      f"min {min(diff):.1f}, max {max(diff):.1f}")
print(f"rescue pairs skipped : {sum(int(r[7]) for r in ok)}")
print(f"genomes with a CHANGED report or gff: {len(changed)}")
for a in changed:
    print(f"  CHANGED {a}")
PYEOF
rm -rf "$WORK/fasta"
echo "done."
