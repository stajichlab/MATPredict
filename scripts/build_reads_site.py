#!/usr/bin/env python3
"""Write the data the raw-reads pages of the site read (web/intro/reads/): one index and one file per species (JSON for the page, TSV to download).

Input: a wave folder of `datasets/reads/sra_survey/pilot/wave1/` (results_all_final.tsv, gates_final.tsv) and the species list (wave1_species.tsv). Standard library only.
    build_reads_site.py --out web/intro/data/reads [--wave datasets/reads/sra_survey/pilot/wave1]
The pages say "confirmed / provisional / withheld"; the wave driver's gate calls them PASS / PROVISIONAL / HOLD.
"""
import argparse, csv, itertools, json, os, re, sys
from collections import Counter

CONF = {"PASS": "confirmed", "PROVISIONAL": "provisional", "PROVISIONAL_HOMOTHALLIC": "provisional (homothallic)"}
FIELDS = ["run", "strain", "study", "layout", "call", "truth", "verdict", "shared_depth", "depth_1", "depth_2", "breadth_1", "breadth_2", "flags", "rescue", "truth_assembly", "clone"]


def rd(p):
    return list(csv.DictReader(open(p, newline="").read().replace("\r", "").splitlines(), delimiter="\t"))


CLONE_ANI, CLONE_ALT, MIN_RUNS_FOR_COUNTS = 0.998, 0.997, 10


def clone_groups(pairs_file, rows):
    """Clone groups of one species from its saved pair identities (clone_pilot/<slug>.pairs.tsv, written by clone_clusters.py --pairs): average linkage at 0.998. Sets row["clone"] to "c<N>" for runs in a group
    of two or more. Returns the summary for the species report. Counts per clone are given only where the number of clones is stable between 0.998 and 0.997 (a plateau: within 10% of the runs) and there are at
    least 10 runs; where it is not, similar runs may be close relatives and not copies, and the report says so."""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "datasets", "reads", "sra_survey", "pilot"))
    import clone_clusters as cc
    prs = [(r["run_a"], r["run_b"], float(r["ani"]) if r["ani"] else None) for r in rd(pairs_file)]
    runs = sorted({x for p in prs for x in p[:2]}); P = {frozenset((a, b)): v for a, b, v in prs if v is not None}
    g1 = sorted(cc.cluster(runs, prs, CLONE_ANI, "average"), key=lambda c: (-len(c), c[0])); g2 = cc.cluster(runs, prs, CLONE_ALT, "average")
    by = {r["run"]: r for r in rows}; groups = []
    for i, c in enumerate(g1, 1):
        if len(c) < 2: continue
        vals = [P[frozenset(q)] for q in itertools.combinations(c, 2) if frozenset(q) in P]
        for x in c:
            if x in by: by[x]["clone"] = f"c{i}"
        groups.append({"id": f"c{i}", "size": len(c), "runs": sorted(c), "calls": dict(Counter(by[x]["call"] for x in c if x in by)), "min_ani": round(min(vals), 4) if vals else None})
    plateau = abs(len(g1) - len(g2)) <= 0.10 * len(runs) and len(runs) >= MIN_RUNS_FOR_COUNTS
    out = {"threshold": CLONE_ANI, "linkage": "average", "runs": len(runs), "clones": len(g1), "clones_at_alt": len(g2), "alt_threshold": CLONE_ALT, "groups": groups, "plateau": plateau, "counts": None}
    if plateau:
        cnt = Counter()
        for c in g1:
            calls = {by[x]["call"] for x in c if x in by and by[x]["call"] in ("MAT1-1", "MAT1-2", "both")}
            if calls: cnt[next(iter(calls)) if len(calls) == 1 else "conflict"] += 1
        out["counts"] = dict(cnt)
    return out


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


SINGLE = ("MAT1-1", "MAT1-2")
CURATED = [
    dict(slug="fola", species="Fusarium oxysporum f. sp. lactucae (Fola)", src="results/2026-10-06_fola_reads_type/v4_concordance.tsv", truth_label="the samtools breadth call (bwa-mem2 per idiomorph, breadth of 90% or more)", tier="own panel",
         study="Fola collection", intro="148 Fola strains with trimmed read pairs, typed from the first 8 million reads with a panel cut from Fola genomes (AT141 and VSP-0980). The comparison call comes from the same reads, so the two methods share any sample mix-up. VSP-0947 is a strain with reads of both idiomorphs.",
         links=[["Analysis note", "analysis/2026-10-06_fola-reads-type.md"]]),
    dict(slug="afum", species="Aspergillus fumigatus (curated collection)", src="results/2026-10-07_afum_reads/comparison.tsv", truth_label="a BLAST of the two idiomorph regions against the assembly of the strain", tier="own panel",
         study="A. fumigatus collection", intro="331 read sets of 304 strains with assemblies plus 28 with reads only; 20 strains have a single-end layout. Truth is a BLAST of the two idiomorph regions against the strain's assembly (assemblies holding both regions, and partial ones, have no single-allele truth). Panel filtered for k-mers that differ between flank alleles (--min-unique-run). The wave 1 report for this species is a separate, smaller sample.",
         links=[["Analysis note", "analysis/2026-10-07_afum-reads-and-assembly.md"]]),
    dict(slug="cocci", species="Coccidioides immitis and C. posadasii", src="datasets/reads/coccidioides/results/cocci2_reads_vs_detect.tsv", truth_label="detect on the strain's assembly (not an independent test)", tier="own panel",
         study="Population-genomics collection", intro="558 strains from the Coccidioides population-genomics collection, typed from the first 8 million reads with one panel from the C. immitis H538.4 (MAT1-1) and RS (MAT1-2) loci. Strains with too little flank depth are not called. Two strains carry corrected species labels. The wave 1 reports for the two species are a separate, smaller sample.",
         links=[["Dataset notes", "datasets/reads/README.md"]]),
]


def curated(c, repo, out):
    rows = rd(os.path.join(repo, c["src"])); R = []
    if c["slug"] == "afum": rows = [r for r in rows if r["reads"] != "no_reads"]      # one strain has no read set
    for r in rows:
        if c["slug"] == "fola": call, truth, st, run, d1, d2, sh, fl = r["kmer_call"], r["truth"], "", r["sample"], r["MAT1-1_depth"], r["MAT1-2_depth"], r["shared_depth"], r["flags"]; strain = r["sample"]
        elif c["slug"] == "afum": call, truth, run, strain, d1, d2, sh, fl = r["reads"], r["truth"], r["read_prefix"], r["strain"], r["m11_depth"], r["m12_depth"], r["shared_depth"], r["flags"]
        else: call, truth, run, strain, d1, d2, sh, fl = r["reads_call"], r["detect_call"], r["strain"], r["strain"], r["MAT1-1_depth"], r["MAT1-2_depth"], r["shared_depth"], "relabelled species" if r["relabelled"] == "yes" else ""
        sp = r["species"] if c["slug"] == "cocci" else ""
        t = truth if truth in SINGLE else ""
        v = "" if not t else "untyped" if call not in ("MAT1-1", "MAT1-2", "both") else "both_vs_one" if call == "both" else "agree" if call == t else "DISAGREE"
        R.append(dict(run=run, strain=strain, study=sp or c["study"], layout=r.get("layout", ""), call=call, truth=truth if truth else "", verdict=v, shared_depth=sh, depth_1=d1, depth_2=d2, breadth_1="", breadth_2="", flags=fl if fl != "-" else "", rescue="", truth_assembly=""))
    cnt = Counter(r["call"] for r in R); single = [r for r in R if r["call"] in SINGLE and r["truth"] in SINGLE]
    s = {"slug": c["slug"], "kind": "curated", "species": c["species"], "tier": c["tier"], "pooled_strains": len(R), "confidence": "curated collection", "why": [], "runs": len(R), "done": len(R),
         "typed": sum(cnt[k] for k in SINGLE), "MAT1-1": cnt["MAT1-1"], "MAT1-2": cnt["MAT1-2"], "both": cnt["both"], "none": cnt["none"], "low_depth": cnt["low_depth"],
         "with_assembly_call": len(single), "agree": sum(r["verdict"] == "agree" for r in single), "disagree": sum(r["verdict"] == "DISAGREE" for r in single), "rescued": 0,
         "intro": c["intro"], "truth_label": c["truth_label"], "links": c["links"], "refused_with_truth": sum(1 for r in R if r["truth"] in SINGLE and r["call"] in ("low_depth", "none")), "both_vs_one": sum(1 for r in R if r["verdict"] == "both_vs_one")}
    json.dump({**s, "gate": "", "rows": R}, open(os.path.join(out, c["slug"] + ".json"), "w"), separators=(",", ":"))
    with open(os.path.join(out, c["slug"] + ".tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, FIELDS, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(R)
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--wave", default="datasets/reads/sra_survey/pilot/wave1"); ap.add_argument("--out", required=True); ap.add_argument("--repo", default="."); ap.add_argument("--clones", default="datasets/reads/sra_survey/pilot/clone_pilot", help="folder with <slug>.pairs.tsv files")
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
        rows = [{k: r.get(k, "") for k in FIELDS} for r in runs]
        pf = os.path.join(a.repo, a.clones, slug + ".pairs.tsv"); cl = clone_groups(pf, rows) if os.path.exists(pf) else None
        if cl: s["clone_groups"] = len(cl["groups"]); s["clone_runs"] = cl["runs"]
        index.append(s)
        json.dump({**s, "gate": g["reasons"], "clones": cl, "rows": rows}, open(os.path.join(a.out, slug + ".json"), "w"), separators=(",", ":"))
        with open(os.path.join(a.out, slug + ".tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, FIELDS, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    allr = [r for r in res if r["state"] == "done"]
    typed = [r for r in allr if r["call"] in ("MAT1-1", "MAT1-2")]; tr = [r for r in typed if r["truth"]]
    tot = {"species": len(index), "runs": len(res), "done": len(allr), "typed": len(typed), "with_assembly_call": len(tr), "agree": sum(r["verdict"] == "agree" for r in tr),
           "disagree": sum(r["verdict"] == "DISAGREE" for r in tr), "confidence": dict(Counter(s["confidence"] for s in index)),
           "calls": dict(Counter(r["call"] for r in allr)), "rescued": sum(1 for r in typed if r.get("rescue"))}
    cur = [curated(c, a.repo, a.out) for c in CURATED if os.path.exists(os.path.join(a.repo, c["src"]))]
    json.dump({"totals": tot, "species": index, "curated": cur}, open(os.path.join(a.out, "index.json"), "w"), separators=(",", ":"))
    print(f"{len(index)} species, {tot['runs']} runs; {tot['agree']}/{tot['with_assembly_call']} agree; {tot['confidence']}")


if __name__ == "__main__":
    main()
