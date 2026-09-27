"""Before/after comparison for the Rhodotorula redPR/redHD records.

Usage: compare.py TRUTH_TSV   (truth = strain allele calls derived from the
unpublished manuscript's tables; kept outside the repository)
Writes per_genome.tsv and prints the summary.
"""
import collections, csv, glob, os, sys
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
RECORD_ASM = {"GCA_921037615.3", "GCA_000988875.2", "GCA_920103745.3", "GCA_024748845.1",
              "GCA_060419215.1", "GCA_002917965.1", "GCA_026119225.1"}  # strains that supplied a record
truth = {}
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    if r["asmid"]:
        truth[r["asmid"]] = r
samples = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}


def load(arm):
    out = {}
    for f in glob.glob(f"{HERE}/run_{arm}/runs/*/detection_report.yaml"):
        asm = f.split("/")[-2]
        r = yaml.safe_load(open(f)) or {}
        calls = r.get("detected") or []
        wall = open(os.path.join(os.path.dirname(f), "wall_seconds")).read().strip() if os.path.exists(
            os.path.join(os.path.dirname(f), "wall_seconds")) else ""
        out[asm] = dict(calls=calls, wall=float(wall) if wall else None, routing=r.get("routing_mode"))
    return out


before, after = load("ac33880"), load("02434ee")
rows = []
for asm in sorted(set(before) | set(after)):
    a = after.get(asm, {"calls": []}); b = before.get(asm, {"calls": []})
    def fam(calls, f):
        return [c for c in calls if c["family"].endswith(":" + f)]
    pr_a = sorted({c["idiomorph"] for c in fam(a["calls"], "redPR")})
    pr_b = sorted({c["idiomorph"] for c in fam(b["calls"], "redPR")})
    hd_a = fam(a["calls"], "redHD")
    pr_contigs = {c["contig"] for c in fam(a["calls"], "redPR")}
    hd_contigs = {c["contig"] for c in hd_a}
    t = truth.get(asm, {})
    s = samples.get(asm, {})
    rows.append(dict(asmid=asm, species=s.get("SPECIES", ""), held_out=asm.split("_")[0] + "_" + asm.split("_")[1] not in RECORD_ASM,
                     has_known_answer=bool(t),
                     PR_agrees_known=("" if not t else str("+".join(pr_a) == t.get("PR", "").replace("/", "+"))),
                     PR_before="+".join(pr_b), PR_after="+".join(pr_a),
                     HD_after=len(hd_a), HD_conf="+".join(sorted({c["confidence"] for c in hd_a})),
                     HD_PR_same_contig=bool(pr_contigs & hd_contigs),
                     wall_before=b.get("wall"), wall_after=a.get("wall")))
with open(f"{HERE}/per_genome.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

n = len(rows)
print(f"genomes {n}; before PR called {sum(1 for r in rows if r['PR_before'])}; after PR called {sum(1 for r in rows if r['PR_after'])}; "
      f"after HD called {sum(1 for r in rows if r['HD_after'])}")
print("HD and PR on the same contig (after):", sum(1 for r in rows if r["HD_PR_same_contig"]))
# known answers
for label, sel in (("all", lambda r: True), ("held-out", lambda r: r["held_out"])):
    k = [r for r in rows if r["has_known_answer"] and sel(r)]
    agree = collections.Counter()
    for r in k:
        t = truth[r["asmid"]]["PR"].replace("/", "+")
        agree[("after", "agree" if r["PR_after"] == t else ("uncalled" if not r["PR_after"] else "disagree"))] += 1
        agree[("before", "agree" if r["PR_before"] == t else ("uncalled" if not r["PR_before"] else "disagree"))] += 1
        agree[("HD", "called" if r["HD_after"] else "uncalled")] += 1
    print(f"known answers ({label}) n={len(k)}:", dict(agree))
    for r in k:
        if r["PR_after"] != truth[r["asmid"]]["PR"].replace("/", "+"):
            print("   ", r["asmid"], r["species"], "known answer is a hybrid" if "/" in truth[r["asmid"]]["PR"] else "", "after", r["PR_after"] or "-", "before", r["PR_before"] or "-")
# species allele split (after)
sp = collections.defaultdict(collections.Counter)
for r in rows:
    sp[r["species"]][r["PR_after"] or "uncalled"] += 1
print("allele split by species (after):")
for k2, v in sorted(sp.items(), key=lambda kv: -sum(kv[1].values())):
    print("   ", k2, dict(v))
wb = sorted(r["wall_before"] for r in rows if r["wall_before"]); wa = sorted(r["wall_after"] for r in rows if r["wall_after"])
if wb and wa:
    print(f"median wall s: before {wb[len(wb)//2]:.0f}, after {wa[len(wa)//2]:.0f}")
