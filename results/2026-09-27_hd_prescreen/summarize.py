import csv, statistics as st, collections
rows = [r for r in csv.DictReader(open("scores.tsv"), delimiter="\t") if r["record_genome"] != "True"]
for r in rows:
    for k in ("score", "hd1", "hd2", "max_blast_bits"): r[k] = float(r[k])
pos = [r for r in rows if r["label"] == "pos"]
nadm = [r for r in rows if r["label"] == "neg" and r["admitted"] == "True"]
nnon = [r for r in rows if r["label"] == "neg" and r["admitted"] != "True"]
def q(v): v = sorted(v); return (round(v[0],1), round(st.median(v),1), round(v[-1],1)) if v else None
print("evaluated (record genomes excluded): pos", len(pos), "neg admitted", len(nadm), "neg non-admitted sample", len(nnon))
print("HMM score min/median/max  pos", q([r["score"] for r in pos]), " neg-admitted", q([r["score"] for r in nadm]), " neg-nonadm", q([r["score"] for r in nnon]))
print("tblastn best bits          pos", q([r["max_blast_bits"] for r in pos]), " neg-admitted", q([r["max_blast_bits"] for r in nadm]))
worst_pos = min(r["score"] for r in pos)
for name, key in (("HMM", "score"), ("blast", "max_blast_bits")):
    wp = min(r[key] for r in pos)
    removed = sum(1 for r in nadm if r[key] < wp)
    print(f"{name}: threshold = worst positive {wp:.1f}; admitted negatives removed {removed}/{len(nadm)} ({100*removed/len(nadm):.0f}%); best admitted negative {max(r[key] for r in nadm):.1f}")
# negatives above worst positive: what are they
above = sorted([r for r in nadm if r["score"] >= worst_pos], key=lambda r: -r["score"])
print("admitted negatives scoring >= worst positive:", len(above))
for r in above[:12]: print("  ", r["species"], r["genome"][:30], r["score"], "blast", r["max_blast_bits"], "capped", r["polish_capped"])
lo = sorted(pos, key=lambda r: r["score"])[:8]
print("lowest positives:"); [print("  ", r["species"], r["genome"][:30], r["score"], "hd1", r["hd1"], "hd2", r["hd2"], "LOSO-excl", r["loso_excluded"]) for r in lo]
# per-genome polish reduction
g_adm = collections.Counter(r["genome"] for r in rows if r["admitted"] == "True")
g_keep = collections.Counter(r["genome"] for r in rows if r["admitted"] == "True" and r["score"] >= worst_pos)
print("admitted clusters per genome median", st.median(g_adm.values()), "-> kept median", st.median([g_keep.get(g,0) for g in g_adm]))
