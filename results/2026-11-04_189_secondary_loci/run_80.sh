set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict; O=$R/results/2026-11-04_189_secondary_loci/run; mkdir -p $O/lists
python3 - <<'P'
import csv
R="/bigdata/stajichlab/jstajich/projects/MATPredict/"
t={r["ASMID"]:r["NCBI_TAXONID"] for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
g=[l.strip() for l in open(R+"results/2026-11-04_189_secondary_loci/g189.txt")]
open(R+"results/2026-11-04_189_secondary_loci/run/lists/g80.tsv","w").write("".join(f"{a}\t{t[a]}\n" for a in g if a in t)); print(sum(a in t for a in g),"of",len(g))
P
W=$R/.claude/worktrees/run-flags
DETECT_ARGS="--preset v0.8.1" sbatch --parsable -p short --time=2:00:00 --mem=24gb --cpus-per-task=8 -J i189 -o $R/logs/i189_%j.log \
 --export=ALL,REPO=$R,SRC=$W/src,MATPREDICT_DB_ROOT=$W/db,CLADE=i189,SAMPLE_LIST=$O/lists/g80.tsv,OUT_DIR=$O/out,GENOME_TIMEOUT=3600 $W/scripts/run_clade_panel.slurm
