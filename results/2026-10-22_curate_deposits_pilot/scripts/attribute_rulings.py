"""Attribute v1 -> v2 verdict changes to the 2026-10-22 curator rulings."""
import csv, json, sys, collections
v1 = {r["protein"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
v2 = json.load(open(sys.argv[2]))
out = sys.argv[3]

def rulings(checks):
    rs = []
    for c in checks:
        if c["status"] != "needs_human":
            continue
        n, d = c["name"], c["detail"]
        if n == "hold": rs.append("1_hold")
        elif n == "other_idiomorph": rs.append("4_two_gene")
        elif n in ("conflict",) or (n == "duplicate" and "same model" in d): rs.append("8_curated_conflict_or_duplicate")
        elif n == "intron" and ("GC-AG" in d or "shorter than" in d): rs.append("7_doubtful_intron")
        elif n == "gene_span": rs.append("7_wide_gene_feature")
        elif n in ("start", "stop") and "partial" in d: rs.append("10_partial")
        elif n == "species": rs.append("13_species_not_binomial")
        else: rs.append("other:" + n)
    return sorted(set(rs))

trans = collections.Counter(); per = collections.Counter(); sole = collections.Counter(); all_nh = collections.Counter()
rows = []
for r in v2:
    p = r["protein"]; b = v1[p]["verdict"]; a = r["verdict"]
    trans[(b, a)] += 1
    rs = rulings(r["checks"])
    for x in rs: all_nh[x] += 1
    moved = b != "NEEDS_HUMAN" and a == "NEEDS_HUMAN"
    new_rs = [x for x in rs if not x.startswith("other:")]
    if moved:
        for x in new_rs: per[x] += 1
        if len(new_rs) == 1: sole[new_rs[0]] += 1
    f = r["facts"]
    rows.append([p, f.get("organism", ""), f.get("species", ""), f.get("order", ""), f.get("gene", ""), b, a,
                 "yes" if moved else "", ",".join(rs)])
with open(out, "w") as fh:
    fh.write("protein\torganism\tspecies\torder\tgene\tverdict_v1\tverdict_v2\tmoved_to_needs_human\tneeds_human_causes\n")
    for row in rows: fh.write("\t".join(row) + "\n")
print("transitions v1 -> v2:")
for k, v in sorted(trans.items()): print(f"  {k[0]:22s} -> {k[1]:22s} {v}")
print("deposits moved to NEEDS_HUMAN, per ruling (a deposit may count under several):")
for k, v in sorted(per.items()): print(f"  {k:34s} {v}  (sole cause: {sole.get(k, 0)})")
print("all NEEDS_HUMAN causes in v2 (deposit counts):")
for k, v in sorted(all_nh.items()): print(f"  {k:34s} {v}")
