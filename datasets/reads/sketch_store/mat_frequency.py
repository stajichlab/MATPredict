#!/usr/bin/env python3
"""MAT idiomorph frequency of a sketch store, per sample and per clone group, overall and by study, country, source and lineage.

    mat_frequency.py --store STORE [--query STORE/query] [--clone 0.9995] [--lineage 0.998] [--min-mean-ani 0.95] --out OUTDIR

Inputs (all written by the store build): QUERY/mat/mat_calls.tsv (sketch and reads-type calls), QUERY/distances/samples_qc.tsv
(mean ANI, outlier flag, clone groups at several thresholds), QUERY/qc/qc.tsv (k-mer QC flags), STORE/exclude.tsv (curator
exclusions, optional) and STORE/manifest.tsv (strain, first_public).

Call per sample: the reads-type call when it is MAT1-1 or MAT1-2, else the sketch call when that is single; 'both' (both idiomorphs,
a mixture or a heterothallic-looking sample) and unresolved calls are counted but not in the frequency.
Samples left out of the frequency, each with its reason in samples.tsv: curator exclusion; mean ANI below --min-mean-ani (another
species: an A. fumigatus store at 0.95 keeps divergent lineages, mean ANI ~0.99, and drops section Fumigati relatives, ~0.83);
a k-mer QC 'mixture?' flag.
Clone-corrected frequency: one vote per clone group (average-linkage groups at --clone ANI), the group's majority call (ties: no
vote), so a study that sequenced one clone 50 times counts once. Lineages are the groups at --lineage ANI.
Fractions are MAT1-1 / (MAT1-1 + MAT1-2) with a Wilson 95% interval and a two-sided exact binomial test against 0.5.
Writes OUTDIR/samples.tsv, OUTDIR/freq_<grouping>.tsv and OUTDIR/summary.txt."""
import argparse, csv, collections, math, os, sys

SINGLE = ("MAT1-1", "MAT1-2")


def rd(path):
    if not path or not os.path.exists(path):
        return {}
    return {r["sample_id"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def wilson(k, n, z=1.96):
    if not n:
        return (float("nan"),) * 2
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def binom_p(k, n):
    """Two-sided exact binomial test of k successes in n against p = 0.5."""
    if not n:
        return float("nan")
    pmf = [math.comb(n, i) / 2 ** n for i in range(n + 1)]
    return min(1.0, sum(x for x in pmf if x <= pmf[k] * (1 + 1e-9)))


def call(m):
    rt, sk = m.get("reads_type_call", ""), m.get("sketch_call", "")
    if rt in SINGLE:
        return rt
    if sk in SINGLE:
        return sk
    return "both" if "both" in (rt, sk) else "unresolved"


def freq_rows(groups):
    """groups: {name: Counter of calls} -> rows with n, MAT1-1, MAT1-2, fraction, CI, p."""
    out = []
    for g, c in sorted(groups.items(), key=lambda x: -(x[1]["MAT1-1"] + x[1]["MAT1-2"])):
        a, b = c["MAT1-1"], c["MAT1-2"]; n = a + b; lo, hi = wilson(a, n)
        out.append(dict(group=g, n_typed=n, MAT1_1=a, MAT1_2=b, both=c["both"], unresolved=c["unresolved"],
                        frac_MAT1_1=f"{a / n:.3f}" if n else "", ci95_lo=f"{lo:.3f}" if n else "", ci95_hi=f"{hi:.3f}" if n else "",
                        p_vs_0_5=f"{binom_p(a, n):.3g}" if n else ""))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--store", required=True); ap.add_argument("--query"); ap.add_argument("--out", required=True)
    ap.add_argument("--clone", default="0.9995"); ap.add_argument("--lineage", default="0.998")
    ap.add_argument("--min-mean-ani", type=float, default=0.95); ap.add_argument("--min-group", type=int, default=10,
                    help="groups with fewer typed samples are pooled as 'other' in the by-group tables")
    a = ap.parse_args(argv)
    q = a.query or os.path.join(a.store, "query")
    mat, sq, qc = rd(f"{q}/mat/mat_calls.tsv"), rd(f"{q}/distances/samples_qc.tsv"), rd(f"{q}/qc/qc.tsv")
    excl, man = rd(os.path.join(a.store, "exclude.tsv")), rd(os.path.join(a.store, "manifest.tsv"))
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for s, m in mat.items():
        d, f = sq.get(s, {}), qc.get(s, {})
        reason = ""
        if s in excl:
            reason = "curator_exclude"
        elif d.get("mean_ani") and float(d["mean_ani"]) < a.min_mean_ani:
            reason = "other_species"
        elif "mixture?" in (f.get("flag") or ""):
            reason = "qc_mixture"
        country = (m.get("country") or "").split(":")[0].strip() or ("local" if m.get("source") == "local" else "unknown")
        rows.append(dict(sample_id=s, source=m.get("source", ""), study=m.get("study_accession", ""), country=country,
                         strain=man.get(s, {}).get("strain", ""), call=call(m), excluded=reason, mean_ani=d.get("mean_ani", ""),
                         clone=d.get(f"clone_{a.clone}", ""), lineage=d.get(f"clone_{a.lineage}", "")))
    with open(f"{a.out}/samples.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    kept = [r for r in rows if not r["excluded"]]
    # one vote per clone group: majority single call; a sample without a clone group (outlier, unsketched) votes alone
    votes = collections.defaultdict(collections.Counter)
    for r in kept:
        votes[r["clone"] or "solo:" + r["sample_id"]][r["call"]] += 1
    clone_rows = []
    for g, c in votes.items():
        a1, a2 = c["MAT1-1"], c["MAT1-2"]
        v = "MAT1-1" if a1 > a2 else "MAT1-2" if a2 > a1 else ("unresolved" if not a1 else "tie")
        first = next(r for r in kept if (r["clone"] or "solo:" + r["sample_id"]) == g)
        clone_rows.append(dict(first, call=v if v != "tie" else "unresolved", n_members=sum(c.values())))
    summary = []

    def table(name, items, key, pool):
        g = collections.defaultdict(collections.Counter)
        for r in items:
            g[key(r)][r["call"]] += 1
        if pool:
            small = [k for k, c in g.items() if c["MAT1-1"] + c["MAT1-2"] < a.min_group]
            for k in small:
                g["other"].update(g.pop(k))
        fr = freq_rows(g)
        with open(f"{a.out}/freq_{name}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, list(fr[0]) if fr else ["group"], delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(fr)
        return fr

    ex = collections.Counter(r["excluded"] or "kept" for r in rows)
    summary.append(f"samples with a MAT row: {len(rows)}; {dict(ex)}")
    for label, items in (("samples", kept), ("clone-corrected", clone_rows)):
        fr = table(f"overall_{label}", items, lambda r: "all", False)[0]
        summary.append(f"{label:16s} MAT1-1 {fr['MAT1_1']}  MAT1-2 {fr['MAT1_2']}  both {fr['both']}  unresolved {fr['unresolved']}  "
                       f"MAT1-1 fraction {fr['frac_MAT1_1']} (95% CI {fr['ci95_lo']}-{fr['ci95_hi']}, p vs 0.5 {fr['p_vs_0_5']})")
    summary.append(f"clone groups at ANI {a.clone}: {len(clone_rows)} (from {len(kept)} samples)")
    for name, key in (("source", lambda r: r["source"]), ("study", lambda r: r["study"]), ("country", lambda r: r["country"]),
                      ("lineage", lambda r: r["lineage"] or "none")):
        for label, items in (("samples", kept), ("clones", clone_rows)):
            fr = table(f"{name}_{label}", items, key, name != "source")
            if label == "clones":
                sig = [f"{r['group']} {r['frac_MAT1_1']} (n={r['n_typed']}, p={r['p_vs_0_5']})" for r in fr if r["p_vs_0_5"] and float(r["p_vs_0_5"]) < 0.01]
                summary.append(f"by {name} (clone-corrected, groups >= {a.min_group} typed): {len(fr)} groups; "
                               f"differ from 0.5 at p<0.01: {'; '.join(sig[:12]) or 'none'}")
    open(f"{a.out}/summary.txt", "w").write("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    sys.exit(main())
