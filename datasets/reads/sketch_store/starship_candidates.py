#!/usr/bin/env python3
"""Per-isolate Starship tables and the tiered bafA / Nebuchadnezzar-h1 candidate list from a store's screen.

    starship_candidates.py --screen STORE/query/screen_starships_bafA/screen.tsv --samples matfreq/samples.tsv --out OUTDIR

--samples is the samples.tsv written by mat_frequency.py (study, source, country, strain, excluded, clone, lineage per sample).
Writes OUTDIR/starship_calls_wide.tsv (every sketched isolate x every target: hash fraction and a P/p/- call),
OUTDIR/starship_summary.tsv (per element: kept isolates present / partial / absent, clone groups present by majority vote) and
OUTDIR/bafA_neb_h1_candidates_v2.tsv (tiers below).

Tiers: strong = both remaining bafA hashes found (fraction 1.0) and Nebuchadnezzar-h1 present; probable = bafA present and
Nebuchadnezzar-h1 present or partial; discordant = one of the two present without the other. 'excluded' carries the sample flag
(other species, k-mer mixture, curator exclusion). Local isolates appear here by name: do not commit this table; the markdown report
uses only its public rows (starship_report_md.py)."""
import argparse, collections, csv, sys

LETTER = {"present": "P", "partial": "p", "absent": "-"}
TIER_ORDER = {"strong": 0, "probable": 1, "discordant": 2}


def rd(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(rows)


def tier_of(bf, bc, nc):
    if bf >= 1.0 and nc == "present":
        return "strong"
    if bc == "present" and nc in ("present", "partial"):
        return "probable"
    if bc == "present" or nc == "present":
        return "discordant"
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--screen", required=True); ap.add_argument("--samples", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    samples = {r["sample_id"]: r for r in rd(a.samples)}
    scr, targets = collections.defaultdict(dict), []
    for r in rd(a.screen):
        scr[r["sample_id"]][r["target"]] = (float(r["fraction"]), r["call"])
        if r["target"] not in targets:
            targets.append(r["target"])
    for need in ("bafA", "Nebuchadnezzar-h1"):
        if need not in targets:
            sys.exit(f"target {need} is not in the screen")
    ids = sorted(scr)
    write(f"{a.out}/starship_calls_wide.tsv", ["sample_id", "strain", "source", "study", "country", "clone", "lineage", "excluded"] + [f"{t}_frac" for t in targets] + [f"{t}_call" for t in targets],
          [[s] + [samples.get(s, {}).get(k, "") for k in ("strain", "source", "study", "country", "clone", "lineage", "excluded")]
           + [f"{scr[s][t][0]:.3f}" for t in targets] + [LETTER[scr[s][t][1]] for t in targets] for s in ids])
    kept = [s for s in ids if s in samples and not samples[s]["excluded"]]
    summ = []
    for t in targets:
        c = collections.Counter(scr[s][t][1] for s in kept)
        g = collections.defaultdict(collections.Counter)
        for s in kept:
            g[samples[s]["clone"] or "solo:" + s][scr[s][t][1]] += 1
        cp = sum(1 for v in g.values() if v["present"] > v["absent"] and v["present"] >= v["partial"])
        summ.append([t, len(kept), c["present"], c["partial"], c["absent"], len(g), cp, f"{cp / len(g):.3f}" if g else ""])
    write(f"{a.out}/starship_summary.tsv", ["element", "kept_isolates", "present", "partial", "absent", "clone_groups", "clones_present", "clone_frac_present"], summ)
    out, counts = [], collections.Counter()
    for s in ids:
        bf, bc = scr[s]["bafA"]; nf, nc = scr[s]["Nebuchadnezzar-h1"]
        t = tier_of(bf, bc, nc)
        if not t:
            continue
        m = samples.get(s, {})
        counts[(t, "excluded" if m.get("excluded") else "kept")] += 1
        out.append([t, s, m.get("strain", ""), m.get("source", ""), m.get("study", ""), m.get("country", ""), f"{bf:.2f}", bc, f"{nf:.3f}", nc, m.get("excluded", ""), m.get("clone", ""), m.get("lineage", "")])
    out.sort(key=lambda x: (TIER_ORDER[x[0]], x[11], x[1]))
    write(f"{a.out}/bafA_neb_h1_candidates_v2.tsv", ["tier", "sample_id", "strain", "source", "study", "country", "bafA_frac", "bafA_call", "neb_h1_frac", "neb_h1_call", "excluded", "clone", "lineage"], out)
    print(sorted(counts.items()))
    print("strong kept clone groups", len({x[11] or x[1] for x in out if x[0] == "strong" and not x[10]}))


if __name__ == "__main__":
    sys.exit(main())
