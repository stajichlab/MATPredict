#!/usr/bin/bash
# HPCC: build the Ogataea sample list (134 genomes) and a DB copy that holds the six candidate MATog records.
set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict
O=$R/results/2026-11-04_matog_candidates/loo; mkdir -p $O/lists $O/rep
python3 - <<'P' > $O/acc.txt
import csv
for r in csv.DictReader(open("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-31_v080_ascomycota/genomes.tsv"), delimiter="\t"):
    if r["genus"] == "Ogataea": print(r["genome"])
P
sed 's#^#*/#;s#$#/detection_report.yaml#' $O/acc.txt > $O/pat.txt
cd $O/rep && zstd -dc $R/results/2026-10-31_v080_ascomycota/reports_all.tar.zst | tar -x --wildcards -T $O/pat.txt
python3 - <<'P'
import glob, yaml
O = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-04_matog_candidates/loo/"
rows = []
for f in sorted(glob.glob(O + "rep/**/detection_report.yaml", recursive=True)):
    d = yaml.safe_load(open(f)); rows.append((f.split("/")[-2], d["run"]["taxid"]))
open(O + "lists/all.tsv", "w").write("".join(f"{g}\t{t}\n" for g, t in rows)); print(len(rows), "genomes")
P
rm -rf $O/db_cand; cp -r $R/.claude/worktrees/run-a2933f7/db $O/db_cand
for d in $R/results/2026-11-04_matog_candidates/records_dirs/*; do cp -r $d $O/db_cand/Ascomycota/Pichiales/; done
ls $O/db_cand/Ascomycota/Pichiales | grep -c MATog
