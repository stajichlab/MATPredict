#!/usr/bin/bash -l
#SBATCH -N 1 -c 1 --mem 2gb --time 00:20:00 -J doth-taxonkit -p short
# Verify the Dothideomycetes taxid and find BFD genomes under it whose CLASS column is blank or differs.
set -euo pipefail
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test/bin:$PATH
D=$1
CSV=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/ncbi_accessions_taxonomy.csv
cd "$D"
echo 147541 | taxonkit lineage -n -r > taxid_147541.txt
cat taxid_147541.txt
taxonkit list --ids 147541 -n -r --indent "  " 2>/dev/null | awk '/order|subclass/' > dothideomycetes_orders_taxonkit.txt || true
python3 - "$CSV" <<'PY'
import csv,sys
r=list(csv.DictReader(open(sys.argv[1])))
open("all_asco_taxids.txt","w").write("\n".join(sorted({x["NCBI_TAXID"] for x in r if x["NCBI_TAXID"]}))+"\n")
PY
taxonkit lineage all_asco_taxids.txt -n -r > all_taxid_lineage.tsv 2>taxonkit.err
grep -c . all_taxid_lineage.tsv
awk -F'\t' '$2 ~ /Dothideomycetes/ {print $1}' all_taxid_lineage.tsv > taxids_under_147541.txt
wc -l taxids_under_147541.txt
