"""Split-locus rule (4174440) vs the same code without it (8d80bed).

For each panel: calls gained, lost, changed; every split_locus call listed with
its core gene, identity, edge distance, flanks and contigs, plus the
classifier verdict. Paralog-looking calls are flagged: a split call whose
core identity is < 97, whose flanks sit > 5 points below the core, or where
the genome already has a call of another family on the flank contigs.
"""
import csv
import os
import sys

import yaml

R = os.path.dirname(os.path.abspath(__file__))
BASE = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
samp = {r["ASMID"]: r for r in csv.DictReader(open(SAMPLES))}

PANELS = [
    ("Mucoromycota", f"{BASE}/2026-09-27_flank_bitscore_implemented/Mucoromycota_8d80bed",
     f"{R}/Mucoromycota_4174440"),
    ("Cryptococcus", f"{R}/Cryptococcus_8d80bed", f"{R}/Cryptococcus_4174440"),
    ("Serinales200", f"{R}/Serinales200_8d80bed", f"{R}/Serinales200_4174440"),
    ("Dothideomycetes", f"{R}/Dothideomycetes_8d80bed", f"{R}/Dothideomycetes_4174440"),
]


def load(d):
    out = {}
    runs = os.path.join(d, "runs")
    if not os.path.isdir(runs):
        return out
    for a in os.listdir(runs):
        f = os.path.join(runs, a, "detection_report.yaml")
        if not os.path.exists(f):
            continue
        r = yaml.safe_load(open(f)) or {}
        out[a] = r.get("detected") or []
    return out


def key(c):
    return (c["family"], c["contig"], c["idiomorph"])


def main():
    lines = []
    for name, before_dir, after_dir in PANELS:
        b, a = load(before_dir), load(after_dir)
        common = sorted(set(b) & set(a))
        gained = lost = changed = 0
        split = []
        for g in common:
            bk = {(c["family"], c["contig"]): c for c in b[g]}
            ak = {(c["family"], c["contig"]): c for c in a[g]}
            for k in set(ak) - set(bk):
                gained += 1
            for k in set(bk) - set(ak):
                lost += 1
            for k in set(ak) & set(bk):
                if (ak[k]["idiomorph"], ak[k]["confidence"], ak[k].get("locus_class")) != \
                   (bk[k]["idiomorph"], bk[k]["confidence"], bk[k].get("locus_class")):
                    changed += 1
            for c in a[g]:
                if c.get("split_locus"):
                    split.append((g, c, len(b[g])))
        genomes_called_b = sum(1 for g in common if b[g])
        genomes_called_a = sum(1 for g in common if a[g])
        lines.append(f"## {name}: genomes {len(common)} (before {len(b)}, after {len(a)}); "
                     f"called {genomes_called_b} -> {genomes_called_a}; calls gained {gained}, "
                     f"lost {lost}, changed {changed}; split_locus calls {len(split)}")
        for g, c, nb in split:
            s = c["split_locus"]
            sp = samp.get(g, {})
            clf = c.get("idiomorph_classifier") or {}
            flag = []
            if s["core_identity"] < 97:
                flag.append("core<97")
            if any(f["identity"] < s["core_identity"] - 5 for f in s["flanks"]):
                flag.append("flank>5below")
            if nb:
                flag.append("genome-had-other-call")
            fl = ", ".join(f"{f['gene']}@{f['contig']}:{f['identity']}" for f in s["flanks"])
            lines.append(
                f"   {g} | {sp.get('SPECIES', '')} {sp.get('STRAIN', '')} | {c['family']} "
                f"{c['idiomorph']} {c['confidence']} | core {s['core_gene']} {s['core_identity']} "
                f"on {s['core_contig']} edge {s['edge_distance']} | flanks {fl} | "
                f"classifier {clf.get('input', clf.get('classifier_input', ''))} "
                f"margin {clf.get('margin', '')} | {'FLAG ' + ','.join(flag) if flag else 'ok'}"
            )
    text = "\n".join(lines)
    open(os.path.join(R, "compare_output.txt"), "w").write(text + "\n")
    print(text)


if __name__ == "__main__":
    sys.exit(main())
