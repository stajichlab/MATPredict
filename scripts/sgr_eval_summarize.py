"""Summarise the strong-single-gene evaluation. Usage: sgr_eval_summarize.py RUN_DIR TRUTH_TSV"""
import csv, glob, json, sys
from collections import Counter
from pathlib import Path

run = Path(sys.argv[1]); truth = {r["genome"]: r["stated"] for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
lists = {p.stem: {l.split("\t")[0] for l in p.read_text().splitlines() if l} for p in (run / "lists").glob("*.tsv")}
rows = {}
for f in glob.glob(str(run / (sys.argv[3] if len(sys.argv) > 3 else "out") / "summary.*.jsonl")):
    for l in open(f):
        r = json.loads(l); rows[r["genome"]] = r
print("genomes done", len(rows), "errors", sum("error" in r for r in rows.values()))
sig = lambda l: (l["contig"], l["start"], l["end"], l["idiomorph"], l["confidence"], l["locus_class"], l["family"])
def added(r):
    off = {sig(l) for l in r["off"]["loci"]}
    return [l for l in r["on"]["loci"] if sig(l) not in off]
def changed(r):  # off loci that are missing or altered in on
    on = {sig(l) for l in r["on"]["loci"]}
    return [l for l in r["off"]["loci"] if sig(l) not in on]
def agree(gid, l):
    t = truth.get(gid)
    if t not in ("MAT1-1", "MAT1-2") or l["family"] != "Ascomycota:MAT": return None
    return l["idiomorph"] == t
out = []
for name in ("named5", "pleo_uncalled", "pleo_called", "control_called", "uncalled_nonpleo", "truth203"):
    ids = [g for g in lists.get(name, ()) if g in rows and "error" not in rows[g]]
    n_off_called = sum(bool(rows[g]["off"]["loci"]) for g in ids)
    adds = {g: added(rows[g]) for g in ids}
    gen_add = [g for g in ids if adds[g]]
    gen_rec = [g for g in gen_add if not rows[g]["off"]["loci"]]
    gen_extra = [g for g in gen_add if rows[g]["off"]["loci"]]
    chg = [g for g in ids if changed(rows[g])]
    print(f"\n== {name}: {len(ids)}/{len(lists.get(name, ()))} done; called off: {n_off_called}")
    print(f"   genomes with loci added: {len(gen_add)} (recovered uncalled {len(gen_rec)}; extra locus in an already-called genome {len(gen_extra)}); loci added {sum(len(v) for v in adds.values())}")
    unflag = [(g, l) for g in ids for l in adds[g] if not l["strong"]]
    print(f"   UNFLAGGED added loci (not made by the rule): {len(unflag)} {[(g[:15], l['idiomorph'], l['confidence']) for g, l in unflag][:5]}")
    print(f"   genomes whose off loci changed/lost under on: {len(chg)} {chg[:5]}")
    n2 = sum(rows[g]['on']['two_idiomorphs'] > rows[g]['off']['two_idiomorphs'] for g in ids)
    print(f"   genomes with a new two_idiomorphs statement: {n2}")
    ag = Counter()
    for g in gen_add:
        for l in adds[g]:
            a = agree(g, l); ag["no_truth" if a is None else ("agree" if a else "DISAGREE")] += 1
    print("   added loci vs truth:", dict(ag))
    for g in gen_add:
        for l in adds[g]:
            s = l["strong"] or {}
            print(f"   +{g[:15]} {l['family']} {l['idiomorph']} {l['confidence']}/{l['locus_class']} model_id={s.get('identity')} cov={s.get('coverage')} hit_id={s.get('cluster_identity')} off_called={bool(rows[g]['off']['loci'])} truth={truth.get(g)}")
    if name == "truth203":
        for tag in ("off", "on"):
            c = Counter()
            for g in ids:
                labs = {l["idiomorph"] for l in rows[g][tag]["loci"] if l["family"] == "Ascomycota:MAT" and l["idiomorph"] in ("MAT1-1", "MAT1-2")}
                t = truth[g]
                c["no_call" if not labs else ("agree" if labs == {t} else ("both" if len(labs) > 1 and t in labs else "DISAGREE"))] += 1
            print(f"   truth203 {tag}:", dict(c))
tm = [(g, r["off_s"], r["on_s"]) for g, r in rows.items() if "off_s" in r]
missing = [g for n in lists for g in lists[n] if g not in rows]
print("\nlisted but not in results:", len(set(missing)), sorted(set(missing))[:5])
print("\nmean off s", sum(t[1] for t in tm) / max(1, len(tm)), "mean on s", sum(t[2] for t in tm) / max(1, len(tm)))
