set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict; O=$R/results/2026-11-04_148_hyphopichia; P=$R/results/2026-10-24_ctg_mtla1_records; mkdir -p $O/lists
cat > $O/build_testdb.py <<'P'
import shutil, sys, yaml
from pathlib import Path
sys.path.insert(0, "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-cassette/src")
R = Path("/bigdata/stajichlab/jstajich/projects/MATPredict"); O = R / "results/2026-11-04_148_hyphopichia"; P = R / "results/2026-10-24_ctg_mtla1_records"
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.cli import _http_transport, build_gff_for_record
from MATPredict.db.build_duckdb import build as build_duckdb
ncbi = NcbiClient(email=None, api_key=None, fetcher=CachedFetcher(cache_dir=P / "ncbi_cache", transport=_http_transport))
t = O / "testdb"
if t.exists(): shutil.rmtree(t)
shutil.copytree(R / ".claude/worktrees/run-a2933f7/db", t)
c = P / "cand_hyph/Ascomycota/717740_makgeolli_MTL_A"
dst = t / "Ascomycota/Serinales" / c.name; shutil.copytree(c, dst)
m = yaml.safe_load((dst / "metadata.yaml").read_text()); m["validation"]["status"] = "accepted"
(dst / "metadata.yaml").write_text(yaml.safe_dump(m, sort_keys=False, width=100))
build_gff_for_record(t, "Ascomycota", "Serinales", c.name, ncbi, None)
print(c.name, sorted(p.name for p in dst.iterdir()), flush=True)
build_duckdb(db_root=t, out_path=t / "matpredict.duckdb"); print("duckdb built")
P
$R/.pixi/envs/default/bin/python $O/build_testdb.py
python3 - <<'P'
import csv, random
R="/bigdata/stajichlab/jstajich/projects/MATPredict/"; O=R+"results/2026-11-04_148_hyphopichia/lists/"
t={r["ASMID"]:r["NCBI_TAXONID"] for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
named="GCA_030557965.1 GCA_030566005.1 GCA_030571075.1 GCA_030578995.1 GCA_030580355.1 GCF_019202705.1 GCA_030444895.1 GCA_030556945.1 GCA_030563265.1 GCA_030567015.1 GCA_030674215.1 GCA_003856795.1".split()
ser=[r["genome"] for r in csv.DictReader(open(R+"results/2026-10-31_v080_ascomycota/genomes.tsv"),delimiter="\t") if r["order"]=="Serinales" and "MTL" in r["families_called"]]
random.seed(148); pick=random.sample(ser,min(380,len(ser)))
allg=[g for g in ser if any(g.startswith(n) for n in named)]+pick
seen=[]; [seen.append(g) for g in allg if g not in seen]
open(O+"s400.tsv","w").write("".join(f"{g}\t{t[g]}\n" for g in seen if g in t)); print(len(seen),"genomes of",len(ser),"Serinales MTL-called")
P
W=$R/.claude/worktrees/run-cassette
for tag in base:$R/.claude/worktrees/run-a2933f7/db cand:$O/testdb; do n=${tag%%:*}; d=${tag#*:}
DETECT_ARGS="--preset v0.8.1" sbatch --parsable -p short --time=4:00:00 --mem=24gb --cpus-per-task=8 -J h148-$n -o $R/logs/h148_${n}_%j.log \
 --export=ALL,REPO=$R,SRC=$W/src,MATPREDICT_DB_ROOT=$d,CLADE=h148-$n,SAMPLE_LIST=$O/lists/s400.tsv,OUT_DIR=$O/$n,GENOME_TIMEOUT=3600 $W/scripts/run_clade_panel.slurm; done
