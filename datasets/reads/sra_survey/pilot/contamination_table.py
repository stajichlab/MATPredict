#!/usr/bin/env python3
"""One table of the runs the contamination screen flagged (GmG on the typing subsample; screen_reads.py), across every committed wave folder.

status mixed = GmG's own mixture flag (a second genus at about 2% of the pairs at 2M spots); genus_mismatch = the genus the run was selected for is
neither GmG's call nor among its top hits (wrong organism or mislabelled); borderline = a second genus with >= 5 hashes and concentration >= 0.7 that
GmG did not flag (screened again at 6M spots). The run is typed either way; these are flags for a curator, not exclusions.

    contamination_table.py PILOT_DIR OUT.tsv
"""
import csv, glob, os, re, sys
PILOT, OUT = sys.argv[1:3]
FLAGGED = ("mixed", "genus_mismatch", "borderline")
rows = []
for w in sorted(glob.glob(os.path.join(PILOT, "wave*"))):
    wave = os.path.basename(w)
    if not re.fullmatch(r"wave\d+[a-z_]*", wave) or not os.path.exists(os.path.join(w, "results_all_final.tsv")): continue
    sp = {g["slug"]: g["species"] for g in csv.DictReader(open(os.path.join(w, "gates_final.tsv"), newline=""), delimiter="\t")} if os.path.exists(os.path.join(w, "gates_final.tsv")) else {}
    for r in csv.DictReader(open(os.path.join(w, "results_all_final.tsv"), newline=""), delimiter="\t"):
        if r.get("screen") in FLAGGED:
            rows.append({"run": r["run"], "species": sp.get(r["group"], r["group"]), "strain": r.get("strain", ""), "study": r.get("study", ""),
                         "screen": r["screen"], "screen_genus": r.get("screen_genus", ""), "second_genus": r.get("screen_second", ""),
                         "second_hashes": r.get("screen_hashes", ""), "second_share": r.get("screen_explained", ""), "deeper_screen": r.get("screen_deeper", ""),
                         "call": r.get("call", ""), "truth": r.get("truth", ""), "wave": wave})
rows.sort(key=lambda x: (x["screen"] != "genus_mismatch", x["species"], x["run"]))
cols = ["run", "species", "strain", "study", "screen", "screen_genus", "second_genus", "second_hashes", "second_share", "deeper_screen", "call", "truth", "wave"]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
from collections import Counter
print(len(rows), "flagged runs;", dict(Counter(r["screen"] for r in rows)), "; species:", len({r["species"] for r in rows}))
for s, n in Counter((r["species"], r["screen"]) for r in rows).most_common(15): print(" ", s, n)
