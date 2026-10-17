"""Sample genome lists for the strong-single-gene evaluation (seed 1). Run on the HPCC."""
import csv, random, sys
from collections import defaultdict
from pathlib import Path

CAMP = Path("/bigdata/stajichlab/jstajich/projects/mat-prof/results/2026-10-03_ascomycota_v060")
TRUTH = Path(sys.argv[1])  # truth_set.tsv
OUT = Path(sys.argv[2]); OUT.mkdir(exist_ok=True, parents=True)
FA = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes")
tax = {}
for f in ("GCA", "GCF"):
    for line in open(f"/bigdata/stajichlab/jstajich/projects/biosample_truth/assembly_summary_{f}.txt"):
        if line.startswith("#"): continue
        c = line.split("\t"); tax[c[0]] = c[5]
rows = list(csv.DictReader(open(CAMP / "genomes.tsv"), delimiter="\t"))
def ok(r):
    acc = "_".join(r["genome"].split("_")[:2])
    return acc in tax and (FA / (r["genome"] + ".fa.gz")).exists()
rows = [r for r in rows if r["order"] and ok(r)]
def line(gid):
    acc = "_".join(gid.split("_")[:2]); return f"{gid}\t{tax[acc]}\t{FA/(gid+'.fa.gz')}"
rng = random.Random(1)
def w(name, ids):
    (OUT / name).write_text("".join(line(g) + "\n" for g in ids)); print(name, len(ids))
# control: called, ten largest orders
cnt = defaultdict(int)
for r in rows: cnt[r["order"]] += 1
top = [o for o, _ in sorted(cnt.items(), key=lambda t: -t[1])[:10]]
print("top orders", top)
ctrl = []
for o in top:
    pool = sorted(r["genome"] for r in rows if r["order"] == o and r["status"] == "called")
    ctrl += rng.sample(pool, min(60, len(pool)))
w("control_called.tsv", ctrl)
pool = sorted(r["genome"] for r in rows if r["status"] == "uncalled" and r["order"] != "Pleosporales")
w("uncalled_nonpleo.tsv", rng.sample(pool, 400))
for st, n in (("uncalled", 24), ("called", 24)):
    pool = sorted(r["genome"] for r in rows if r["order"] == "Pleosporales" and r["status"] == st)
    w(f"pleo_{st}.tsv", rng.sample(pool, n))
named = ["GCA_003497985.1", "GCA_003503065.1", "GCA_022559985.1", "GCA_003503165.1", "GCA_035078975.1"]
allg = {"_".join(r["genome"].split("_")[:2]): r["genome"] for r in rows}
w("named5.tsv", [allg[a] for a in named if a in allg])
t = [r for r in csv.DictReader(open(TRUTH), delimiter="\t") if r["stated"] in ("MAT1-1", "MAT1-2")]
ids = [r["genome"] for r in t if (FA / (r["genome"] + ".fa.gz")).exists() and "_".join(r["genome"].split("_")[:2]) in tax]
w("truth203.tsv", ids)
