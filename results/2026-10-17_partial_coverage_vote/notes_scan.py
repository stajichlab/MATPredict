import sys, glob, os, yaml, csv
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
runs, variant, out = sys.argv[1:4]
w = csv.writer(open(out, "w"), delimiter="\t", lineterminator="\n")
w.writerow(["genome", "contig", "start", "label_on", "confidence", "decision", "winner_gene", "winner_cov", "winner_bits", "other_gene", "other_cov", "other_bits", "edge", "profile_scores", "profile_margin", "reason"])
n = 0
for p in sorted(glob.glob(f"{runs}/{variant}/*/detection_report.yaml")):
    doc = yaml.load(open(p), Loader=L)
    for c in doc.get("detected") or []:
        v = c.get("idiomorph_vote_note")
        if v:
            n += 1
            pr = v.get("profile") or {}
            w.writerow([os.path.basename(os.path.dirname(p)), c["contig"], c["start"], c["idiomorph"], c["confidence"], v["decision"],
                        v["winner_gene"], v["winner_coverage"], v["winner_bitscore"], v["other_gene"], v["other_coverage"], v["other_bitscore"],
                        v.get("winner_edge_distance", ""), pr.get("scores", ""), pr.get("margin_bits", ""), v["reason"]])
print(n, "notes")
