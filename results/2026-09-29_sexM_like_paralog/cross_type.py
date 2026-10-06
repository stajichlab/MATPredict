"""Is a called core gene present in genomes of the OPPOSITE mating type?

For each called core protein (the winning sexM for Minus calls, sexP for Plus calls)
tblastn it against up to K confident opposite-type genomes of the same genus
(confident = only that idiomorph called, at high confidence, locus_class mat_locus).
A true idiomorph gene should have no near-identical copy there; a paralog present in
both mating types will. Writes cross_type.tsv.
"""
import collections, csv, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

GEN = "/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/genomes"
SCR = sys.argv[1]
K = 4
rows = list(csv.DictReader(open("models.tsv"), delimiter="\t"))
by = collections.defaultdict(list)
for r in rows:
    by[r["org"]].append(r)

# per-genome idiomorph summary
genome = {}
for org, rs in by.items():
    calls = {}
    for r in rs:
        calls[r["call"]] = (r["idiomorph"], r["confidence"], r["locus_class"])
    idios = {c[0] for c in calls.values()}
    strong = {c[0] for c in calls.values() if c[1] == "high" and c[2] == "mat_locus"}
    genome[org] = (idios, strong)

def genus(org):
    return org.split("_")[0]

conf_only = collections.defaultdict(lambda: collections.defaultdict(list))
for org, (idios, strong) in genome.items():
    if len(idios) == 1 and idios == strong and list(idios)[0] in ("Plus", "Minus"):
        conf_only[genus(org)][list(idios)[0]].append(org)

# winning core model per call
win = {}
core = collections.defaultdict(list)
for r in rows:
    if r["gene"] in ("sexM", "sexP"):
        core[(r["org"], r["call"])].append(r)
for k, v in core.items():
    idio = v[0]["idiomorph"]
    want = {"Minus": "sexM", "Plus": "sexP"}.get(idio)
    c = [x for x in v if x["gene"] == want]
    if c:
        win[k] = max(c, key=lambda x: float(x["bitscore"] or 0))

os.makedirs(SCR, exist_ok=True)

def db(org):
    p = os.path.join(SCR, org)
    if not os.path.exists(p + ".nsq") and not os.path.exists(p + ".00.nsq"):
        subprocess.run(["makeblastdb", "-in", f"{GEN}/{org}.sorted.fasta", "-dbtype", "nucl",
                        "-out", p], check=True, capture_output=True)
    return p

jobs = []
for (org, call), w in win.items():
    opp = "Plus" if w["idiomorph"] == "Minus" else "Minus"
    targets = [g for g in sorted(conf_only[genus(org)][opp]) if g != org][:K]
    for t in targets:
        jobs.append((org, call, w, t))

needed = sorted({j[3] for j in jobs})
with ThreadPoolExecutor(4) as ex:
    list(ex.map(db, needed))

def run(job):
    org, call, w, t = job
    q = os.path.join(SCR, f"q_{abs(hash((org, call)))}.faa")
    if not os.path.exists(q):
        open(q, "w").write(f">q\n{w['protein']}\n")
    out = subprocess.run(["tblastn", "-query", q, "-db", os.path.join(SCR, t), "-outfmt",
                          "6 pident length qlen evalue bitscore sseqid sstart send",
                          "-max_target_seqs", "5", "-evalue", "1e-5", "-seg", "no"],
                         capture_output=True, text=True).stdout.strip().splitlines()
    best = None
    for l in out:
        f = l.split("\t")
        pid, ln, ql, bs = float(f[0]), int(f[1]), int(f[2]), float(f[4])
        cov = ln / ql
        if best is None or bs > best[3]:
            best = (pid, cov, f[5], bs)
    return dict(org=org, call=call, idiomorph=w["idiomorph"], confidence=w["confidence"],
                locus_class=w["locus_class"], genes_found=w["genes_found"],
                score_minus=w["score_minus"], score_plus=w["score_plus"], core_gene=w["gene"],
                core_len=len(w["protein"]), opp_genome=t,
                best_pid=best[0] if best else 0, best_cov=round(best[1], 2) if best else 0,
                best_bits=best[3] if best else 0, best_contig=best[2] if best else "")

with ThreadPoolExecutor(4) as ex:
    res = list(ex.map(run, jobs))
cols = list(res[0].keys())
with open("cross_type.tsv", "w", newline="") as fo:
    wr = csv.DictWriter(fo, fieldnames=cols, delimiter="\t")
    wr.writeheader()
    wr.writerows(res)
print(len(res), "comparisons;", len(needed), "opposite-type genomes")
