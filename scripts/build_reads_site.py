#!/usr/bin/env python3
"""Write the data the raw-reads pages of the site read (web/intro/reads/): one index and one file per species (JSON for the page, TSV to download).

Input: a wave folder of `datasets/reads/sra_survey/pilot/wave1/` (results_all_final.tsv, gates_final.tsv) and the species list (wave1_species.tsv). Standard library only.
    build_reads_site.py --out web/intro/data/reads [--wave datasets/reads/sra_survey/pilot/wave1]
The pages say "confirmed / provisional / withheld"; the wave driver's gate calls them PASS / PROVISIONAL / HOLD.
"""
import argparse, csv, json, os, re
from collections import Counter

CONF = {"PASS": "confirmed", "PROVISIONAL": "provisional", "PROVISIONAL_HOMOTHALLIC": "provisional (homothallic)"}
FIELDS = ["run", "strain", "study", "layout", "call", "truth", "verdict", "shared_depth", "depth_1", "depth_2", "breadth_1", "breadth_2", "flags", "rescue", "truth_assembly"]


def rd(p):
    return list(csv.DictReader(open(p, newline="").read().replace("\r", "").splitlines(), delimiter="\t"))


def reason(r):
    """Plain wording for one gate reason token, e.g. typed_0.67 -> 'only 67% of runs got a call'."""
    m = re.match(r"([a-z_0-9]+?)_([0-9.]+)(?:of([0-9]+))?$", r)
    if not m: return r.replace("_", " ")
    k, v, n = m.groups()
    if k == "agree": return f"only {v} of {n} runs with an assembly call agree with it"
    if k in ("typed", "provisional_typed"): return f"only {round(float(v) * 100)}% of runs got a single-allele call"
    if k == "shared_depth": return f"flank depth too low ({v}x)"
    if k == "both": return f"{round(float(v) * 100)}% of runs show both idiomorphs"
    if k == "too_few_runs": return f"only {v} runs available"
    return r.replace("_", " ")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--wave", default="datasets/reads/sra_survey/pilot/wave1"); ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    res = rd(os.path.join(a.wave, "results_all_final.tsv")); gates = {r["slug"]: r for r in rd(os.path.join(a.wave, "gates_final.tsv"))}
    sp = {r["slug"]: r for r in rd(os.path.join(os.path.dirname(a.wave.rstrip("/")), "wave1_species.tsv"))}
    by = {}
    for r in res:      # the verdict is recomputed here: the results table keeps the pre-rescue verdicts of runs typed by a later pass
        r["verdict"] = "" if not r["truth"] else "untyped" if r["state"] != "done" or r["call"] not in ("MAT1-1", "MAT1-2", "both") else "both_vs_one" if r["call"] == "both" else "agree" if r["call"] == r["truth"] else "DISAGREE"
        by.setdefault(r["group"], []).append(r)
    index = []
    for slug, g in sorted(gates.items(), key=lambda kv: kv[1]["species"]):
        runs = by.get(slug, []); done = [r for r in runs if r["state"] == "done"]; c = Counter(r["call"] for r in done)
        typed = [r for r in done if r["call"] in ("MAT1-1", "MAT1-2")]
        truth = [r for r in typed if r["truth"]]
        d = g["decision"].split(":")[0]; why = list(dict.fromkeys(reason(x) for x in g["decision"].split(":", 1)[1].split("+"))) if ":" in g["decision"] else []
        s = {"slug": slug, "species": g["species"], "tier": sp.get(slug, {}).get("reference_tier", ""), "pooled_strains": int(sp.get(slug, {}).get("strains_pooled") or 0),
             "confidence": "withheld" if d == "HOLD" else CONF.get(d, d), "why": why, "runs": len(runs), "done": len(done), "typed": len(typed),
             "MAT1-1": c["MAT1-1"], "MAT1-2": c["MAT1-2"], "both": c["both"], "none": c["none"], "low_depth": c["low_depth"],
             "with_assembly_call": len(truth), "agree": sum(r["verdict"] == "agree" for r in truth), "disagree": sum(r["verdict"] == "DISAGREE" for r in truth),
             "rescued": sum(1 for r in typed if r.get("rescue"))}
        index.append(s)
        rows = [{k: r.get(k, "") for k in FIELDS} for r in runs]
        json.dump({**s, "gate": g["reasons"], "rows": rows}, open(os.path.join(a.out, slug + ".json"), "w"), separators=(",", ":"))
        with open(os.path.join(a.out, slug + ".tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, FIELDS, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    allr = [r for r in res if r["state"] == "done"]
    typed = [r for r in allr if r["call"] in ("MAT1-1", "MAT1-2")]; tr = [r for r in typed if r["truth"]]
    tot = {"species": len(index), "runs": len(res), "done": len(allr), "typed": len(typed), "with_assembly_call": len(tr), "agree": sum(r["verdict"] == "agree" for r in tr),
           "disagree": sum(r["verdict"] == "DISAGREE" for r in tr), "confidence": dict(Counter(s["confidence"] for s in index)),
           "calls": dict(Counter(r["call"] for r in allr)), "rescued": sum(1 for r in typed if r.get("rescue"))}
    json.dump({"totals": tot, "species": index}, open(os.path.join(a.out, "index.json"), "w"), separators=(",", ":"))
    print(f"{len(index)} species, {tot['runs']} runs; {tot['agree']}/{tot['with_assembly_call']} agree; {tot['confidence']}")


if __name__ == "__main__":
    main()
