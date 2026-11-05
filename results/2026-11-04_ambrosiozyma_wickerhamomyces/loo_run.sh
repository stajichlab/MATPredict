#!/usr/bin/bash
# HPCC: v0.8.1 on the 18 Ambrosiozyma and 22 Wickerhamomyces genomes with the base DB and with the two candidate families and records (<= 24 CPUs).
set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict; H=$R/results/2026-11-04_ambrosiozyma_wickerhamomyces; O=$H/loo; mkdir -p $O/lists
python3 - <<'P'
import csv
R="/bigdata/stajichlab/jstajich/projects/MATPredict/"; H=R+"results/2026-11-04_ambrosiozyma_wickerhamomyces/loo/lists/"
t={r["ASMID"]:r["NCBI_TAXONID"] for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
g=[r["genome"] for r in csv.DictReader(open(R+"results/2026-10-31_v080_ascomycota/genomes.tsv"),delimiter="\t") if r["genus"] in ("Ambrosiozyma","Wickerhamomyces")]
src={"GCA_003707675.3":"self0","GCA_054949465.1":"self1"}
rows=[(a,t[a]) for a in g if a in t]
open(H+"all.tsv","w").write("".join(f"{a}\t{x}\n" for a,x in rows))
open(H+"others.tsv","w").write("".join(f"{a}\t{x}\n" for a,x in rows if a[:15] not in src))
for a,x in rows:
    if a[:15] in src: open(H+src[a[:15]]+".tsv","w").write(f"{a}\t{x}\n")
print(len(rows),"genomes")
P
rm -rf $O/db_cand; cp -r $R/.claude/worktrees/run-a2933f7/db $O/db_cand
cat $H/families_proposal.yml >> $O/db_cand/Ascomycota/order.yml
for r in 43984_nrrl-y-7523_MATam_MATa-MATalpha:Pichiales 640698_tbrc-bcc18952_MATwe_MATa-MATalpha:Phaffomycetales; do rid=${r%%:*}; ord=${r##*:}; mkdir -p $O/db_cand/Ascomycota/$ord/$rid; cp $H/records/$rid.yaml $O/db_cand/Ascomycota/$ord/$rid/metadata.yaml; cp $H/records/$rid.proteins.faa $O/db_cand/Ascomycota/$ord/$rid/proteins.faa; done
W=$R/.claude/worktrees/run-cassette
sub() { DETECT_ARGS="--preset v0.8.1 $5" sbatch --parsable -p short --time=2:00:00 --mem=24gb --cpus-per-task=$4 -J aw-$1 -o $R/logs/aw_$1_%j.log \
  --export=ALL,REPO=$R,SRC=$W/src,MATPREDICT_DB_ROOT=$2,CLADE=aw-$1,SAMPLE_LIST=$3,OUT_DIR=$O/$1,GENOME_TIMEOUT=3600 $W/scripts/run_clade_panel.slurm; }
sub base $R/.claude/worktrees/run-a2933f7/db $O/lists/all.tsv 8 ""
sub cand $O/db_cand $O/lists/others.tsv 8 ""
sub self0 $O/db_cand $O/lists/self0.tsv 2 "--exclude-records 43984_nrrl-y-7523_MATam_MATa-MATalpha"
sub self1 $O/db_cand $O/lists/self1.tsv 2 "--exclude-records 640698_tbrc-bcc18952_MATwe_MATa-MATalpha"
