#!/usr/bin/env bash
# Download the ENA FASTQ mates of ONE run listed in selected.tsv (one run per Slurm array task).
#
# Usage (Slurm array; task N = Nth data row of selected.tsv, header excluded):
#   N=$(( $(wc -l < selected.tsv) - 1 ))
#   sbatch --array=1-${N}%8 -p short -c 1 --mem=2G -t 04:00:00 -J fus_dl \
#          --wrap "bash $PWD/download.sh selected.tsv /path/to/reads/fusarium"
#   # re-submit only unfinished runs: just run the same command again (idempotent);
#   # or a subset: --array=1-100%8   (rows 1-100; selected.tsv is ordered group/species)
# Usage (single task / testing):
#   SLURM_ARRAY_TASK_ID=3 bash download.sh selected.tsv /path/to/outdir
#   bash download.sh selected.tsv /path/to/outdir SRR28734943      # by run accession
#
# Arguments: 1 = TSV (selected.tsv or manifest.tsv; needs columns run_accession, fastq_ftp, md5)
#            2 = target directory (created); files land in <outdir>/<run>/<run>_{1,2}.fastq.gz
#            3 = optional run accession (overrides SLURM_ARRAY_TASK_ID)
# Behaviour: skips a run when <outdir>/<run>/.done exists; resumes partial files (curl -C -);
#   verifies md5 of each file against ENA before renaming .part -> final; retries 5x with backoff.
#   Exit 0 = done/skipped, 1 = failed (re-run to resume).  Keep concurrency modest (%8) -- ENA FTP is shared.
set -euo pipefail

tsv=${1:?usage: download.sh selected.tsv outdir [run_accession]}
outdir=${2:?usage: download.sh selected.tsv outdir [run_accession]}
want=${3:-}

[[ -r $tsv ]] || { echo "cannot read $tsv" >&2; exit 2; }

# pick the row
if [[ -n $want ]]; then
    row=$(awk -F'\t' -v r="$want" 'NR>1 && $1==r' "$tsv")
else
    : "${SLURM_ARRAY_TASK_ID:?set SLURM_ARRAY_TASK_ID or pass a run accession as arg 3}"
    row=$(awk -F'\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "$tsv")
fi
[[ -n $row ]] || { echo "no matching row in $tsv" >&2; exit 2; }

# locate columns by header name so column order may change
col() { awk -F'\t' -v c="$1" 'NR==1{for(i=1;i<=NF;i++) if($i==c){print i; exit}}' "$tsv"; }
ci_run=$(col run_accession); ci_ftp=$(col fastq_ftp); ci_md5=$(col md5)
[[ -n $ci_run && -n $ci_ftp && -n $ci_md5 ]] || { echo "missing run_accession/fastq_ftp/md5 column" >&2; exit 2; }

run=$(cut -f"$ci_run" <<<"$row")
ftp=$(cut -f"$ci_ftp" <<<"$row")
md5s=$(cut -f"$ci_md5" <<<"$row")
[[ -n $ftp && -n $md5s ]] || { echo "$run: no paired FASTQ listed in ENA (see notes column); skipping" >&2; exit 0; }

dest=$outdir/$run
mkdir -p "$dest"
if [[ -e $dest/.done ]]; then echo "$run: already complete"; exit 0; fi

if command -v md5sum >/dev/null 2>&1; then
    md5of() { md5sum "$1" | cut -d' ' -f1; }
else
    md5of() { md5 -q "$1"; }
fi

fetch() {  # url dest.part
    local url=$1 part=$2
    if command -v curl >/dev/null 2>&1; then
        curl -fsS -L --retry 3 --connect-timeout 30 -C - -o "$part" "$url"
    else
        wget -q -c -O "$part" "$url"
    fi
}

IFS=';' read -r -a urls <<<"$ftp"
IFS=';' read -r -a sums <<<"$md5s"
[[ ${#urls[@]} -eq ${#sums[@]} ]] || { echo "$run: url/md5 count mismatch" >&2; exit 1; }

for i in "${!urls[@]}"; do
    u=https://${urls[$i]#*://}          # ENA lists host/path without scheme
    f=$dest/$(basename "$u")
    if [[ -s $f && $(md5of "$f") == "${sums[$i]}" ]]; then continue; fi
    ok=0
    for attempt in 1 2 3 4 5; do
        fetch "$u" "$f.part" || true
        if [[ -s $f.part && $(md5of "$f.part") == "${sums[$i]}" ]]; then
            mv "$f.part" "$f"; ok=1; break
        fi
        echo "$run: attempt $attempt failed for $(basename "$f") (size/md5 mismatch or transfer error)" >&2
        # a complete-but-wrong file will never fix itself by resuming: start over once it looks full
        [[ $attempt -ge 2 ]] && rm -f "$f.part"
        sleep $((attempt * 20))
    done
    [[ $ok -eq 1 ]] || { echo "$run: giving up on $(basename "$f")" >&2; exit 1; }
done

date -u +%FT%TZ > "$dest/.done"
echo "$run: done"
