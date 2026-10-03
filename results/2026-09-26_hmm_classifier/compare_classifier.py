"""Mucoromycota calls: 4a24ffb (model-pair decision) vs 7277e10 (HMM classifier).

usage: compare_classifier.py OLD_RUNS NEW_RUNS > compare_output.txt
Writes label_changes.tsv (every call whose idiomorph changed, with the
classifier's scores) and classifier_calls.tsv (every new call's classifier
verdict). Calls are paired within a genome by contig overlap.
"""
import collections, csv, os, sys
import yaml
from yaml import CSafeLoader as L

OLD, NEW = sys.argv[1], sys.argv[2]
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}


def calls(runs, g):
    p = f"{runs}/{g}/detection_report.yaml"
    if not os.path.exists(p):
        return None
    return (yaml.load(open(p), Loader=L) or {}).get("detected") or []


def ov(a, b):
    return a["contig"] == b["contig"] and not (a["end"] < b["start"] or a["start"] > b["end"])


def clf(b):
    c = b.get("idiomorph_classifier") or {}
    s = c.get("scores") or {}
    return c.get("idiomorph", ""), s.get("Plus", ""), s.get("Minus", ""), c.get("margin", ""), c.get("proteins_scored", "")


genomes = sorted(set(os.listdir(OLD)) & set(os.listdir(NEW)))
tot = collections.Counter()
fam = collections.defaultdict(collections.Counter)
changes, allnew = [], []
for g in genomes:
    o, n = calls(OLD, g), calls(NEW, g)
    if o is None or n is None:
        tot["missing_report"] += 1
        continue
    f = meta.get(g, {}).get("FAMILY", "?")
    sp = meta.get(g, {}).get("SPECIES", "")
    st = meta.get(g, {}).get("STRAIN", "")
    used = set()
    for a in o:
        m = next((i for i, b in enumerate(n) if i not in used and ov(a, b)), None)
        if m is None:
            k = "lost"; b = None
        else:
            used.add(m); b = n[m]
            k = "unchanged" if a["idiomorph"] == b["idiomorph"] else f"{a['idiomorph']}->{b['idiomorph']}"
        tot[k] += 1; fam[f][k] += 1
        if k != "unchanged":
            v = clf(b) if b else ("", "", "", "", "")
            changes.append(dict(genome=g, family=f, species=sp, strain=st, change=k, contig=a["contig"],
                                start=a["start"], old_genes=",".join(a.get("genes_found") or []),
                                new_genes=",".join((b or {}).get("genes_found") or []),
                                clf_verdict=v[0], clf_plus=v[1], clf_minus=v[2], clf_margin=v[3],
                                proteins=v[4], new_conf=(b or {}).get("confidence", "")))
    for i, b in enumerate(n):
        v = clf(b)
        allnew.append(dict(genome=g, family=f, species=sp, strain=st, contig=b["contig"], start=b["start"],
                           idiomorph=b["idiomorph"], confidence=b["confidence"], clf_verdict=v[0],
                           clf_plus=v[1], clf_minus=v[2], clf_margin=v[3], proteins=v[4],
                           genes=",".join(b.get("genes_found") or [])))
        if i not in used:
            tot["new"] += 1; fam[f]["new"] += 1
            changes.append(dict(genome=g, family=f, species=sp, strain=st, change="new", contig=b["contig"],
                                start=b["start"], old_genes="", new_genes=",".join(b.get("genes_found") or []),
                                clf_verdict=v[0], clf_plus=v[1], clf_minus=v[2], clf_margin=v[3],
                                proteins=v[4], new_conf=b["confidence"]))

for name, rows in (("label_changes.tsv", changes), ("classifier_calls.tsv", allnew)):
    if rows:
        with open(name, "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

print(f"genomes compared: {len(genomes)}")
print("overall:", dict(tot))
basis = collections.Counter("classifier" if r["clf_verdict"] else "fallback (no modelled core protein)" for r in allnew)
print("new-code calls by basis:", dict(basis))
print("classifier verdicts:", dict(collections.Counter(r["clf_verdict"] for r in allnew if r["clf_verdict"])))
margins = sorted(float(r["clf_margin"]) for r in allnew if r["clf_verdict"] not in ("", "undetermined"))
if margins:
    print(f"decisive margins: n={len(margins)} min={margins[0]} median={margins[len(margins)//2]}")
print("\nper family (only families with a change):")
for f, c in sorted(fam.items()):
    if any(k != "unchanged" for k in c):
        print(f"  {f:28s} {dict(c)}")
print("\nevery change:")
for r in changes:
    print(f"  {r['change']:22s} {r['species']} {r['strain']} ({r['genome']}) {r['contig']}:{r['start']} "
          f"clf={r['clf_verdict']} Plus={r['clf_plus']} Minus={r['clf_minus']} margin={r['clf_margin']} "
          f"genes {r['old_genes']} -> {r['new_genes']}")
print("\nundetermined:")
for r in allnew:
    if r["clf_verdict"] == "undetermined":
        print(f"  {r['species']} {r['strain']} ({r['genome']}) {r['contig']}:{r['start']} "
              f"Plus={r['clf_plus']} Minus={r['clf_minus']} margin={r['clf_margin']} genes {r['genes']}")
print("\nUmbelopsis and Lichtheimiaceae/Syncephalastraceae calls:")
for r in allnew:
    if r["species"].startswith("Umbelopsis") or r["family"] in ("Lichtheimiaceae", "Syncephalastraceae"):
        print(f"  {r['family']:18s} {r['species']} {r['strain']} {r['idiomorph']:12s} {r['confidence']:6s} "
              f"clf={r['clf_verdict'] or 'none'} Plus={r['clf_plus']} Minus={r['clf_minus']} margin={r['clf_margin']}")
