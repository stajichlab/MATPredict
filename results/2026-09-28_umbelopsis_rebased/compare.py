"""Compare the rebased curation-umbelopsis arms (guard on / guard off) with the
076afe4-equivalent baseline (results/2026-09-28_next_fixes/Mucoromycota_a9fb69c,
no Umbelopsis/S. racemosum records).

Writes calls_<arm>.tsv, guard_effect.tsv, family_table.tsv and prints a summary.
"""
import csv, glob, os, sys
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
ARMS = {
    "base": os.path.join(RES, "2026-09-28_next_fixes/Mucoromycota_a9fb69c"),
    "guard": os.path.join(HERE, "Mucoromycota_667252a"),
    "noguard": os.path.join(HERE, "Mucoromycota_umb-noguard"),
}
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
samp = {r["ASMID"]: r for r in csv.DictReader(open(SAMPLES))}
WATCH = {
    "GCA_900175165.2": "R. pusillus FCH_5_7", "GCA_900079185.1": "A. glauca",
    "GCA_977110945.1": "U. ramanniana gzUmbRama1", "GCF_025399195.1": "U. ramanniana AG (Umbra1)",
}


def load(arm_dir):
    calls, supp = {}, {}
    for f in glob.glob(os.path.join(arm_dir, "runs", "*", "detection_report.yaml")):
        asm = f.split(os.sep)[-2]
        try:
            r = yaml.safe_load(open(f)) or {}
        except Exception:
            continue
        rows = []
        for x in r.get("detected") or []:
            clf = x.get("idiomorph_classifier") or {}
            rows.append(dict(contig=x["contig"], start=x["start"], end=x["end"],
                             idiomorph=x.get("idiomorph"), confidence=x.get("confidence"),
                             locus_class=x.get("locus_class"),
                             clf_input=clf.get("input") or clf.get("classifier_input"),
                             margin=clf.get("margin"), genes=",".join(x.get("genes_found") or [])))
        calls[asm] = rows
        supp[asm] = [(s.get("contig"), s.get("start"), s.get("end"), s.get("withheld_reason"),
                      s.get("idiomorph")) for s in (r.get("suppressed_loci") or [])]
    return calls, supp


data = {a: load(d) for a, d in ARMS.items() if os.path.isdir(d)}
for a, (calls, _) in data.items():
    with open(os.path.join(HERE, f"calls_{a}.tsv"), "w") as fo:
        fo.write("asmid\tspecies\tfamily\tcontig\tstart\tend\tidiomorph\tconfidence\tlocus_class\tclf_input\tmargin\tgenes\n")
        for asm, rows in sorted(calls.items()):
            s = samp.get(asm, {})
            for x in rows:
                fo.write("\t".join(str(v) for v in [asm, s.get("SPECIES", ""), s.get("FAMILY", ""),
                         x["contig"], x["start"], x["end"], x["idiomorph"], x["confidence"],
                         x["locus_class"], x["clf_input"], x["margin"], x["genes"]]) + "\n")


def key(x):
    return (x["contig"], x["start"], x["end"])


def summ(a, b):
    ca, cb = data[a][0], data[b][0]
    common = set(ca) & set(cb)
    gained = lost = changed = 0
    ga = sum(1 for g in common if ca[g]); gb = sum(1 for g in common if cb[g])
    detail = []
    for g in sorted(common):
        A = {key(x): x for x in ca[g]}; B = {key(x): x for x in cb[g]}
        for k in B.keys() - A.keys():
            gained += 1; detail.append(("gained", g, B[k]))
        for k in A.keys() - B.keys():
            lost += 1; detail.append(("lost", g, A[k]))
        for k in A.keys() & B.keys():
            if (A[k]["idiomorph"], A[k]["confidence"]) != (B[k]["idiomorph"], B[k]["confidence"]):
                changed += 1; detail.append(("changed", g, (A[k], B[k])))
    return dict(n=len(common), called_a=ga, called_b=gb, gained=gained, lost=lost, changed=changed), detail


out = []
for a, b in [("noguard", "guard"), ("base", "guard"), ("base", "noguard")]:
    if a in data and b in data:
        s, det = summ(a, b)
        out.append(f"{a} -> {b}: genomes {s['n']}, called {s['called_a']} -> {s['called_b']}, "
                   f"calls gained {s['gained']}, lost {s['lost']}, changed {s['changed']}")
        with open(os.path.join(HERE, f"diff_{a}_vs_{b}.tsv"), "w") as fo:
            for kind, g, x in det:
                sp = samp.get(g, {}).get("SPECIES", "")
                fo.write(f"{kind}\t{g}\t{sp}\t{x}\n")

fam_rows = []
for g in sorted(set().union(*[set(d[0]) for d in data.values()])):
    s = samp.get(g, {})
    fam = s.get("FAMILY", "")
    if fam not in ("Umbelopsidaceae", "Syncephalastraceae") and g.split("_")[0] + "_" + g.split("_")[1] not in WATCH \
            and not any(g.startswith(w) for w in WATCH) and "stolonifer" not in s.get("SPECIES", "") \
            and "griseocyanus" not in s.get("SPECIES", ""):
        continue
    row = [g, s.get("SPECIES", ""), s.get("STRAIN", ""), fam]
    for a in ("base", "noguard", "guard"):
        rows = data.get(a, ({}, {}))[0].get(g)
        row.append("NA" if rows is None else ";".join(
            f"{x['idiomorph']}/{x['confidence']}/{x['clf_input'] or '-'}/{x['margin']}" for x in rows) or "none")
    fam_rows.append(row)
with open(os.path.join(HERE, "family_table.tsv"), "w") as fo:
    fo.write("asmid\tspecies\tstrain\tfamily\tbase_076afe4\tnoguard\tguard\n")
    for r in fam_rows:
        fo.write("\t".join(r) + "\n")
print("\n".join(out))
print(f"family/watch table rows: {len(fam_rows)}")
