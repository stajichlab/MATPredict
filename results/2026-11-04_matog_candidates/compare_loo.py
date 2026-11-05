"""Compare the MATog calls of the 134 Ogataea genomes without (base) and with (cand + self) the six candidate records. Run on the HPCC.
python3 compare_loo.py LOO_DIR OUT_TSV"""
import csv, glob, os, sys, yaml, collections
O, OUT = sys.argv[1], sys.argv[2]
SP = {r["genome"]: r["species"] for r in csv.DictReader(open("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-31_v080_ascomycota/genomes.tsv"), delimiter="\t")}

def calls(run):
    """{genome: [(idiomorph, confidence, top idiomorph score, n_genes, contig, start, end)]} for MATog calls of one run folder."""
    out = {}
    for f in glob.glob(f"{O}/{run}/runs/*/detection_report.yaml"):
        g = f.split("/")[-2]; d = yaml.safe_load(open(f)); rows = []
        for l in d.get("detected", []):
            if "MATog" not in str(l.get("family", "")): continue
            sc = (l.get("idiomorph_candidates") or [{}])[0].get("score")
            rows.append((l.get("idiomorph"), l.get("confidence"), round(sc, 1) if sc is not None else None, l.get("polished_genes"), l.get("contig"), l.get("start"), l.get("end")))
        out[g] = (sorted(rows, key=lambda r: (r[0] or "", r[4], r[5])), d.get("fallback_search"), d["run"]["database"]["records"])
    return out

base = calls("base")
new = {}
for run in ["cand"] + [f"self{i}" for i in range(6)]:
    new.update(calls(run))
rows, cat = [], collections.Counter()
for g in sorted(base):
    b = base[g][0]; n = new.get(g, ([], None, None))[0]
    bi = sorted({r[0] for r in b}); ni = sorted({r[0] for r in n})
    kind = ("two" if len(bi) == 2 else "one" if len(bi) == 1 else "none")
    if g not in new: change = "NOT_RUN"
    elif [(r[0], r[1], r[4], r[5], r[6]) for r in b] == [(r[0], r[1], r[4], r[5], r[6]) for r in n]: change = "same_calls"
    elif len(ni) > len(bi): change = "gained_idiomorph"
    elif len(ni) < len(bi): change = "LOST_idiomorph"
    else: change = "calls_differ"
    cat[(kind, change)] += 1
    rows.append([g, SP.get(g, "?"), kind, ",".join(bi), ",".join(ni), change,
                 ";".join(f"{r[0]}:{r[1]}:{r[2]}:{r[3]}g:{r[4]}:{r[5]}-{r[6]}" for r in b), ";".join(f"{r[0]}:{r[1]}:{r[2]}:{r[3]}g:{r[4]}:{r[5]}-{r[6]}" for r in n),
                 "yes" if base[g][1] else "", "yes" if new.get(g, (0, None))[1] else "", new.get(g, (0, 0, ""))[2]])
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genome", "species", "base_cassettes", "base_idiomorphs", "new_idiomorphs", "change", "base_calls", "new_calls", "base_fallback", "new_fallback", "db_records"]); w.writerows(rows)
for k, v in sorted(cat.items()): print(k, v)
