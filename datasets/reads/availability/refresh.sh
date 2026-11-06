#!/bin/bash
# Refresh the reads-availability snapshot on the HPCC: re-inventory the lab's local read data with PopSeq's tool, then join PopSeq's
# newest kingdom run table to what MATPredict has typed. Writes datasets/reads/availability/snapshots/<run table date>/ in this checkout,
# with CHANGES.md against the newest earlier snapshot. Nothing is downloaded and nothing outside WORK and the snapshot folder is written.
#
#   sbatch datasets/reads/availability/refresh.sh          # from the root of a MATPredict checkout on the HPCC
#   bash   datasets/reads/availability/refresh.sh          # same, interactively (about 2 minutes)
#
# Then review SUMMARY.md and CHANGES.md and commit the snapshot (untyped_runs.tsv.gz is not committed; it is regenerated from the inputs).
#SBATCH -p stajichlab -c 1 --mem 8G -t 2:00:00 -J reads_availability -o reads_availability.%j.log
set -euo pipefail
POPSEQ=${POPSEQ:-/bigdata/stajichlab/shared/projects/Fungi/PopSeq}
PY=${PY:-/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python}   # needs PyYAML for --db
REPO=${REPO:-$(pwd)}
HERE=$REPO/datasets/reads/availability
SNAP=$HERE/snapshots
WORK=${WORK:-/bigdata/stajichlab/jstajich/projects/MATPredict/results/$(date +%F)_reads_availability}
[ -f "$HERE/reads_availability.py" ] || { echo "run from the root of a MATPredict checkout (REPO=$REPO)" >&2; exit 2; }

# 1. Local data. PopSeq's inventory is read-only on Population_Genomics; one subproject at a time so one failure (a dangling symlink,
#    2026-11-06) does not stop the rest. Failures are listed in WORK/popseq_inventory/failed.txt.
INV=$WORK/popseq_inventory
mkdir -p "$INV"
: > "$INV/failed.txt"
"$POPSEQ/.venv/bin/python" "$POPSEQ/scripts/popseq_inventory.py" --list | awk 'NF{print $1}' > "$INV/subprojects.txt"
while read -r s; do
    if ! "$POPSEQ/.venv/bin/python" "$POPSEQ/scripts/popseq_inventory.py" --subproject "$s" --out-root "$INV" > "$INV/log.$(echo "$s" | tr / _).txt" 2>&1; then
        echo "$s" >> "$INV/failed.txt"
    fi
done < "$INV/subprojects.txt"
echo "local inventory: $(wc -l < "$INV/subprojects.txt") subprojects, $(wc -l < "$INV/failed.txt") failed"
# a subproject PopSeq lists that local_subprojects.tsv does not map is reported, not guessed
for d in "$INV"/panels/*/; do
    s=$(basename "$d")
    grep -q "^$s	" "$HERE/local_subprojects.tsv" || echo "NOTE: subproject $s is not in local_subprojects.tsv (not counted)"
done

# 2. Public data: PopSeq's newest kingdom run table. The previous snapshot is the newest one with a different date.
RUNS=$(ls "$POPSEQ"/cache/fungi_runs_ena_*.tsv.gz | sort | tail -1)
DATE=$(basename "$RUNS" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}')
PREV=$(ls -d "$SNAP"/*/ 2>/dev/null | sed 's:/$::' | grep -v "/$DATE\$" | sort | tail -1 || true)
echo "PopSeq run table: $RUNS; previous snapshot: ${PREV:-none}"

# 3. Join.
"$PY" "$HERE/reads_availability.py" --popseq "$POPSEQ" --runs "$RUNS" --local-inventory "$INV" --db "$REPO/db" --out "$SNAP" ${PREV:+--previous "$PREV"}
cp "$INV/failed.txt" "$SNAP/$DATE/local_inventory_failed.txt"
echo "done: $SNAP/$DATE (review SUMMARY.md${PREV:+ and CHANGES.md})"
