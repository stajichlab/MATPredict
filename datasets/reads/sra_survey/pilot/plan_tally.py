#!/usr/bin/env python3
"""Why were species held at planning? Tally the plan folders of a wave by status, and for a species that has a panel lineages table by what its assemblies carry (both idiomorphs, one, or none).

    plan_tally.py --wave WAVE_DIR
"""
import argparse, collections, csv, os


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--wave", required=True); a = ap.parse_args()
    W = os.path.join(a.wave, "plan"); cat = collections.defaultdict(list)
    for s in sorted(os.listdir(W)):
        sf = os.path.join(W, s, "status"); st = open(sf).read().strip() if os.path.exists(sf) else "planned"
        if st == "planned": cat["planned"].append(s); continue
        pd = os.path.join(W, s, "panel"); lt = [f for f in os.listdir(pd) if f.endswith("_lineages.tsv")] if os.path.isdir(pd) else []
        if not lt: cat[st + " | no lineages table (no loci of the species' idiomorph names)"].append(s); continue
        rows = list(csv.DictReader(open(os.path.join(pd, lt[0])), delimiter="\t")); names = [c for c in rows[0] if c not in ("lineage", "loci", "kept") and not c.startswith("rep_")] if rows else []
        n = [sum(int(r[c]) for r in rows) for c in names]
        cat[f"{st} | " + ("no loci" if not any(n) else "only " + names[n.index(max(n))] if min(n or [0]) == 0 else "both idiomorphs: " + " ".join(names))].append(s)
    for k, v in sorted(cat.items(), key=lambda kv: -len(kv[1])): print(len(v), k, "->", ", ".join(v[:6]))


if __name__ == "__main__":
    main()
