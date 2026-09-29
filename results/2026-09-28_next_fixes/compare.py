#!/usr/bin/env python3
"""Compare a scan against a baseline scan: every call change, with the cause.

Usage: compare.py BASE_DIR NEW_DIR [FAMILY_PREFIX]
Each dir holds runs/<ASMID>/detection_report.yaml (or a reports.tar.zst).
Calls are matched per genome by family + contig + overlap. For each lost call
it looks for a suppressed_loci entry at the same place in NEW and prints its
withheld_reason (e.g. mat_gene_gate). Confidence and label changes are listed.
"""
import io
import subprocess
import sys
import tarfile
from collections import Counter
from pathlib import Path

import yaml


def load(d):
    d = Path(d)
    reps = {}
    for p in d.glob("runs/*/detection_report.yaml"):
        reps[p.parent.name] = yaml.safe_load(p.read_text())
    if not reps:
        tz = next(d.glob("reports*.tar.zst"), None)
        if tz:
            data = subprocess.run(["zstd", "-dc", str(tz)], capture_output=True).stdout
            with tarfile.open(fileobj=io.BytesIO(data)) as t:
                for m in t.getmembers():
                    if m.name.endswith("detection_report.yaml"):
                        reps[Path(m.name).parent.name] = yaml.safe_load(t.extractfile(m).read())
    return reps


def same(a, b):
    return (a["family"] == b["family"] and a["contig"] == b["contig"]
            and a["start"] <= b["end"] and b["start"] <= a["end"])


def main():
    base, new = load(sys.argv[1]), load(sys.argv[2])
    pref = sys.argv[3] if len(sys.argv) > 3 else ""
    tally = Counter()
    lines = []
    called_base = called_new = 0
    for g in sorted(set(base) & set(new)):
        b = [x for x in (base[g].get("detected") or []) if x["family"].startswith(pref)]
        n = [x for x in (new[g].get("detected") or []) if x["family"].startswith(pref)]
        sup = [x for x in (new[g].get("suppressed_loci") or []) if x["family"].startswith(pref)]
        called_base += bool(b)
        called_new += bool(n)
        for x in b:
            m = next((y for y in n if same(x, y)), None)
            if m is None:
                s = next((y for y in sup if same(x, y)), None)
                why = s.get("withheld_reason") if s else "not reported"
                extra = ""
                if s and s.get("withheld_reason") == "mat_gene_gate":
                    extra = (f" score={s.get('best_score')} input={s.get('classifier_input')}"
                             f" flanks={s.get('supporting_flanks')}")
                tally[f"lost:{why}"] += 1
                lines.append(f"LOST  {g} {x['contig']}:{x['start']} {x['idiomorph']}/"
                             f"{x['confidence']} -> {why}{extra}")
                continue
            if m["confidence"] != x["confidence"]:
                tally[f"conf:{x['confidence']}->{m['confidence']}"] += 1
                lines.append(f"CONF  {g} {x['contig']}:{x['start']} {x['idiomorph']} "
                             f"{x['confidence']}->{m['confidence']}")
            if m["idiomorph"] != x["idiomorph"]:
                tally["label"] += 1
                lines.append(f"LABEL {g} {x['contig']}:{x['start']} {x['idiomorph']}->{m['idiomorph']}")
            if (m.get("verification") or {}).get("status") != (x.get("verification") or {}).get("status"):
                tally["verification"] += 1
        for y in n:
            if not any(same(x, y) for x in b):
                tally["gained"] += 1
                lines.append(f"GAIN  {g} {y['contig']}:{y['start']} {y['idiomorph']}/{y['confidence']}")
    print(f"genomes in both: {len(set(base) & set(new))}; genomes called {called_base} -> {called_new}")
    print("tally:", dict(sorted(tally.items())))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
