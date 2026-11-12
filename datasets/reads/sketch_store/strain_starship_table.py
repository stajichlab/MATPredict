#!/usr/bin/env python3
"""Starship table for a named set of strains (for example the eAF strains, PRJNA991533).

One row per isolate whose strain name starts with --prefix: clone group, lineage, the call for every screened Starship target (P present, p partial, - absent),
the all-alleles bafA hit count from the baf family screen and the bafA class (cluster-linked = AF293-type target present, standalone divergent = only the
all-alleles screen sees bafA, none). Writes calls.tsv, fractions.tsv, summary.tsv and table.md.

Local (Population_Genomics) isolates are counted in the summary but never named in any output.

    strain_starship_table.py --calls starship_calls_wide.tsv --baf-screen STORE/query/screen_baf_family/screen.tsv --prefix eAF \
        --set Osiris-h3=0.7 Janus-h4=0.7 --title "eAF strains (PRJNA991533)" --out eAF_table

--set raises the present cut for an element (as in starship_dynamics.py): a call is P at fraction >= cut, p between 0.2 and the cut, - below 0.2; it only tightens.
"""
import argparse
import collections
import csv
import os

BAF_CUT = 2


def load_baf(path):
    hits = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["target"] == "bafA_all":
            hits[r["sample_id"]] = int(r["hits"])
    return hits


def baf_class(r, hits):
    if r["bafA_call"] == "P":
        return "cluster_linked"
    return "standalone_divergent" if hits.get(r["sample_id"], 0) >= BAF_CUT else "none"


def call(r, e, tuned):
    c = r[e + "_call"]
    if e in tuned and r[e + "_frac"] != "":
        f = float(r[e + "_frac"])
        if c in ("P", "p") or f >= 0.2:
            return "P" if f >= tuned[e] else "p" if f >= 0.2 else "-"
        return "-"
    return c


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True); ap.add_argument("--baf-screen", required=True); ap.add_argument("--prefix", required=True)
    ap.add_argument("--set", nargs="*", default=[]); ap.add_argument("--title", default=None); ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    tuned = {k: float(v) for k, v in (x.split("=") for x in a.set)}
    hits = load_baf(a.baf_screen)
    rows = [r for r in csv.DictReader(open(a.calls), delimiter="\t") if r["strain"].startswith(a.prefix) and not r["excluded"]]
    els = [c[:-5] for c in rows[0] if c.endswith("_call")] if rows else []
    els = sorted(els, key=lambda e: (e != "bafA", e.split("-")[0], e))
    local = [r for r in rows if r["source"] == "local"]
    shown = sorted([r for r in rows if r["source"] != "local"],
                   key=lambda r: (int(r["lineage"]) if r["lineage"].isdigit() else 9999, int(r["clone"]) if r["clone"].isdigit() else 99999, r["strain"]))
    os.makedirs(a.out, exist_ok=True)
    cls = {r["sample_id"]: baf_class(r, hits) for r in rows}
    cs = {r["sample_id"]: {e: call(r, e, tuned) for e in els} for r in rows}
    head = ["strain", "sample_id", "clone", "lineage", "n_present", "present_elements", "bafA_all_hits", "bafA_class"]
    with open(f"{a.out}/calls.tsv", "w", newline="") as f, open(f"{a.out}/fractions.tsv", "w", newline="") as g:
        w = csv.writer(f, delimiter="\t", lineterminator="\n"); w2 = csv.writer(g, delimiter="\t", lineterminator="\n")
        w.writerow(head + els); w2.writerow(["strain", "sample_id", "clone", "lineage"] + els)
        for r in shown:
            c = cs[r["sample_id"]]; pres = [e for e in els if c[e] == "P"]
            w.writerow([r["strain"], r["sample_id"], r["clone"], r["lineage"], len(pres), ";".join(pres), hits.get(r["sample_id"], 0), cls[r["sample_id"]]] + [c[e] for e in els])
            w2.writerow([r["strain"], r["sample_id"], r["clone"], r["lineage"]] + [r[e + "_frac"] for e in els])
    clones = collections.defaultdict(list)
    for r in rows:
        clones[r["clone"] or "solo:" + r["sample_id"]].append(r)
    summ = []
    for e in els:
        n = collections.Counter(cs[r["sample_id"]][e] for r in rows)
        summ.append([e, n["P"], n["p"], n["-"], sum(1 for v in clones.values() if sum(cs[r["sample_id"]][e] == "P" for r in v) * 2 > len(v)), len(clones)])
    kinds = ("cluster_linked", "standalone_divergent", "none")
    bsum = {k: (sum(cls[r["sample_id"]] == k for r in rows),
                sum(1 for v in clones.values() if sum(cls[r["sample_id"]] == k for r in v) * 2 > len(v))) for k in kinds}
    with open(f"{a.out}/summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["element", "isolates_present", "isolates_partial", "isolates_absent", "clone_groups_present", "clone_groups_total"]); w.writerows(summ)
        for k in kinds:
            w.writerow([f"bafA_all:{k}", bsum[k][0], "", len(rows) - bsum[k][0], bsum[k][1], len(clones)])
    title = a.title or a.prefix
    with open(f"{a.out}/table.md", "w") as f:
        f.write(f"# {title}: Starship screen\n\n{len(rows)} isolates ({len(local)} local, counted but not listed), {len(clones)} clone groups (ANI 0.9995). "
                "Calls from the solo-hash rescreen: P present, p partial (an upper bound), - absent."
                + (" Cuts raised: " + ", ".join(f"{k} {v}" for k, v in tuned.items()) + "." if tuned else "") + "\n\n")
        f.write("## bafA, all alleles\n\nThe `bafA` column below is the AF293-type allele (in the hrmA cluster). The all-alleles screen also finds divergent copies that sit outside that cluster:\n\n"
                "| bafA class | Isolates | Clone groups |\n|---|---:|---:|\n")
        names = dict(cluster_linked="cluster-linked (AF293-type target present)", standalone_divergent="standalone divergent (all-alleles screen only)", none="no bafA")
        for k in kinds:
            f.write(f"| {names[k]} | {bsum[k][0]} | {bsum[k][1]} of {len(clones)} |\n")
        f.write("\n## Per element\n\n| Element | Present | Partial | Absent | Clone groups present |\n|---|---:|---:|---:|---:|\n")
        for e, p, pp, ab, cp, ct in summ:
            f.write(f"| {e} | {p} | {pp} | {ab} | {cp} of {ct} |\n")
        f.write("\n## Per strain\n\n| Strain | Lineage | Clone | bafA class | n | Present elements |\n|---|---|---|---|---:|---|\n")
        for r in shown:
            pres = [e for e in els if cs[r["sample_id"]][e] == "P"]
            f.write(f"| {r['strain']} | {r['lineage']} | {r['clone']} | {cls[r['sample_id']].replace('_', ' ')} | {len(pres)} | {', '.join(pres)} |\n")
    print(f"{len(rows)} isolates ({len(local)} local), {len(clones)} clone groups; bafA classes " + ", ".join(f"{k} {bsum[k][0]}" for k in kinds))


if __name__ == "__main__":
    main()
