#!/usr/bin/env python3
"""Choose the runs of one species for a typing test: one run per strain, natural-isolate studies only, runs of a few Gb, with a `detect` truth where one exists.

    select_runs.py --species "Pyricularia oryzae" --group pyricularia --panel-dir DIR --panel-prefix PFX --n-truth 16 --n-free 4 --out runs.tsv \
        --runs ENA_RUNS.tsv --meta ENA_META.tsv --studies ../studies.tsv --assembly-summary asm_summary.txt --genomes genomes.tsv --loci loci.tsv [--include SRR..,SRR..]
Rules: the study is not flagged in studies.tsv (mutant, evolution, cross, not_wgs); the run has at least --min-gb Gb and at most --max-gb Gb (the one nearest --target-gb
per sample is chosen: a few Gb is plenty, a 128 Gb run only costs download time); at most --per-study runs from one study; truth = the single `detect` idiomorph of a BFD assembly
with the same BioSample (balanced between the idiomorphs); --include adds named runs.
"""
import argparse, collections, csv, random, re


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for k in ("species", "group", "panel-dir", "panel-prefix", "out", "runs", "meta", "studies", "assembly-summary", "genomes", "loci"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--n-truth", type=int, default=16); ap.add_argument("--n-free", type=int, default=4); ap.add_argument("--per-study", type=int, default=4)
    ap.add_argument("--min-gb", type=float, default=0.6); ap.add_argument("--max-gb", type=float, default=12); ap.add_argument("--target-gb", type=float, default=3)
    ap.add_argument("--include", default=""); ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    asm = {}
    for l in open(a.assembly_summary):
        if l.startswith("#"): continue
        p = l.rstrip("\n").split("\t")
        if len(p) > 3: asm[p[0]] = p[2]
    g = {r["genome"]: r for r in csv.DictReader(open(a.genomes), delimiter="\t") if r["species"].startswith(a.species)}
    truth = collections.defaultdict(set); gen_of = {}
    for r in csv.DictReader(open(a.loci), delimiter="\t"):
        if r["genome"] in g and r["locus_class"] in ("mat_locus", "partial_locus"):
            bs = asm.get(re.match(r"(GC[AF]_\d+\.\d+)", r["genome"]).group(1))
            if bs: truth[bs].add(r["idiomorph"]); gen_of[bs] = r["genome"]
    flagged = {r["study"] for r in csv.DictReader(open(a.studies), delimiter="\t") if r["species"] == a.species and r["excluded_as"]}
    meta = {m["run_accession"]: m for m in csv.DictReader(open(a.meta), delimiter="\t")}
    inc = {x for x in a.include.split(",") if x}
    best = {}
    for r in csv.DictReader(open(a.runs), delimiter="\t"):
        if not r["scientific_name"].startswith(a.species) or not r["study_accession"]: continue
        if r["study_accession"] in flagged and r["run_accession"] not in inc: continue
        try: gb = int(r["base_count"] or 0) / 1e9
        except ValueError: continue
        if not (a.min_gb <= gb <= a.max_gb) and r["run_accession"] not in inc: continue
        k = r["sample_accession"]; score = (r["run_accession"] not in inc, r["library_layout"] != "PAIRED", abs(gb - a.target_gb))
        if k not in best or score < best[k][0]: best[k] = (score, r)
    runs = [v[1] for v in best.values()]
    rnd = random.Random(a.seed); chosen = []; per = collections.Counter()
    def take(pool, n):
        pool = sorted(pool, key=lambda r: r["run_accession"]); rnd.shuffle(pool); got = 0
        for r in pool:
            if got >= n: break
            if r in chosen or per[r["study_accession"]] >= a.per_study: continue
            chosen.append(r); per[r["study_accession"]] += 1; got += 1
        return got
    for r in runs:
        if r["run_accession"] in inc and r not in chosen: chosen.append(r)
    one = lambda r: truth.get(r["sample_accession"]) if len(truth.get(r["sample_accession"], ())) == 1 else None
    nt = max(0, a.n_truth - sum(1 for r in chosen if one(r)))
    for idio, n in (("MAT1-1", nt // 2), ("MAT1-2", nt - nt // 2)):
        got = take([r for r in runs if one(r) and next(iter(one(r))) == idio], n)
    take([r for r in runs if r["sample_accession"] not in truth], a.n_free + (a.n_truth - sum(1 for r in chosen if one(r))) + 0)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run_accession", "group", "study_accession", "strain", "library_layout", "read_count", "base_count", "truth_detect", "truth_assembly", "panel_dir", "panel_prefix"])
        for r in chosen:
            m = meta.get(r["run_accession"], {}); t = truth.get(r["sample_accession"])
            w.writerow([r["run_accession"], a.group, r["study_accession"], (m.get("strain") or m.get("isolate") or m.get("sample_alias") or "")[:30], r["library_layout"], r["read_count"], r["base_count"],
                        ",".join(sorted(t)) if t else "", gen_of.get(r["sample_accession"], ""), a.panel_dir, a.panel_prefix])
    c = collections.Counter(",".join(sorted(truth.get(r["sample_accession"], []))) or "none" for r in chosen)
    print(a.species, len(chosen), "runs; truth:", dict(c), "; studies", len({r["study_accession"] for r in chosen}))


if __name__ == "__main__":
    main()
