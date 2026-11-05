set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict; O=$R/results/2026-11-05_timing; mkdir -p $O
python3 - <<'P'
import csv, random
R="/bigdata/stajichlab/jstajich/projects/MATPredict/"; O=R+"results/2026-11-05_timing/"
t={r["ASMID"]:r["NCBI_TAXONID"] for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
g=list(csv.DictReader(open(R+"results/2026-10-31_v080_ascomycota/genomes.tsv"),delimiter="\t"))
random.seed(5)
ns=[x["genome"] for x in g if x["routing"]=="not_searched"]; lin=[x["genome"] for x in g if x["routing"]=="lineage" and x["order"] in ("Eurotiales","Hypocreales","Saccharomycetales")]; pf=[x["genome"] for x in g if x["routing"]=="phylum_fallback"]
pick=random.sample(ns,6)+random.sample(lin,10)+random.sample(pf,4)
open(O+"t20.tsv","w").write("".join(f"{a}\t{t[a]}\n" for a in pick if a in t)); print(len(pick))
P
sub() { n=$1; w=$2; p=$3; DETECT_ARGS="--preset $p" sbatch --parsable -p short --time=2:00:00 --mem=16gb --cpus-per-task=2 -J tm-$n -o $R/logs/tm_${n}_%j.log \
 --export=ALL,THREADS=1,REPO=$R,SRC=$R/.claude/worktrees/$w/src,MATPREDICT_DB_ROOT=$R/.claude/worktrees/$w/db,CLADE=tm-$n,SAMPLE_LIST=$O/t20.tsv,OUT_DIR=$O/$n,GENOME_TIMEOUT=3600 $R/.claude/worktrees/run-a2933f7/scripts/run_clade_panel.slurm; }
rm -rf $O/v070 $O/v080; sub v070 run-v070 v0.7; sub v080 run-a2933f7 v0.8
