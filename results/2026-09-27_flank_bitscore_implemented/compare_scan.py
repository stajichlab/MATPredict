"""Compare two detection runs call by call.

Usage: compare_scan.py OLD_RUNS_DIR NEW_RUNS_DIR OUT_TSV
A call is keyed by (genome, family, contig); if a genome has several calls on
one contig they are matched in start order. Reports new, lost, and changed
idiomorph / confidence / locus_class.
"""
import collections, csv, glob, os, sys
import yaml

old_dir, new_dir, out = sys.argv[1:4]


def load(d):
    calls = {}
    for p in glob.glob(os.path.join(d, "*", "detection_report.yaml")):
        g = os.path.basename(os.path.dirname(p))
        r = yaml.safe_load(open(p)) or {}
        per = collections.defaultdict(list)
        for x in r.get("detected") or []:
            per[(g, x["family"], x["contig"])].append(x)
        for k, xs in per.items():
            for i, x in enumerate(sorted(xs, key=lambda y: y["start"])):
                calls[k + (i,)] = x
    return calls


old, new = load(old_dir), load(new_dir)
rows, tally = [], collections.Counter()
for k in sorted(set(old) | set(new)):
    a, b = old.get(k), new.get(k)
    if a and not b:
        kind = "lost"
    elif b and not a:
        kind = "new"
    else:
        diffs = [f for f in ("idiomorph", "confidence", "locus_class") if a.get(f) != b.get(f)]
        kind = "changed:" + "+".join(diffs) if diffs else "same"
    tally[kind] += 1
    if kind != "same":
        rows.append(dict(genome=k[0], family=k[1], contig=k[2], kind=kind,
                         old_idiomorph=(a or {}).get("idiomorph"), new_idiomorph=(b or {}).get("idiomorph"),
                         old_conf=(a or {}).get("confidence"), new_conf=(b or {}).get("confidence"),
                         old_class=(a or {}).get("locus_class"), new_class=(b or {}).get("locus_class"),
                         new_ignored=";".join((b or {}).get("confidence_ignored_genes") or []),
                         genes=";".join((b or a).get("genes_found") or [])))
with open(out, "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]) if rows else ["genome"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"old {len(old)} calls, new {len(new)} calls")
for k, v in sorted(tally.items()):
    print(f"  {k}: {v}")
