"""Per-genome idiomorph set of the 2026-10-05 with-records run (results/2026-10-05_dothideomycetes_full/reports_all.tar.zst,
extracted by the caller into a directory)."""
import csv, glob, os, sys
import yaml
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
src, out = sys.argv[1], sys.argv[2]
rows = []
for p in glob.glob(src + "/*/runs/*/detection_report.yaml"):
    g = os.path.basename(os.path.dirname(p))
    d = yaml.load(open(p), Loader=L)
    det = d.get("detected") or []
    rows.append((g, len(det), "+".join(sorted({str(x["idiomorph"]) for x in det})),
                 ";".join(f"{x['idiomorph']}:{x['confidence']}:{x['locus_class']}" for x in det)))
with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["genome", "n_loci", "idiomorphs", "loci"])
    w.writerows(rows)
print(len(rows), "reports")
