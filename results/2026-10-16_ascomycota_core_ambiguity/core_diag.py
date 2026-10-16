import csv, glob, sys, collections, statistics, yaml, os
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
rows = []
for path in glob.glob(sys.argv[1] + "/**/detection_report.yaml", recursive=True):
    genome = os.path.basename(os.path.dirname(path))
    doc = yaml.load(open(path), Loader=L)
    for c in doc.get("detected") or []:
        if c.get("family") != "Ascomycota:MAT":
            continue
        ev = {}
        for g in c.get("gene_evidence") or []:
            if g["gene"] in ("MAT1-1-1", "MAT1-2-1") and not g.get("lost_to"):
                if g["gene"] not in ev or (g.get("bitscore") or 0) > (ev[g["gene"]].get("bitscore") or 0):
                    ev[g["gene"]] = g
        rows.append({"genome": genome, "contig": c["contig"], "start": c["start"], "idiomorph": c["idiomorph"],
                     "locus_class": c["locus_class"], "confidence": c["confidence"],
                     "a": ev.get("MAT1-1-1"), "h": ev.get("MAT1-2-1")})
print("loci", len(rows))
def f(g, k):
    return None if g is None else g.get(k)
# thresholds from single-gene loci (one core gene only): these are the unambiguous positives
single_a = [r["a"] for r in rows if r["a"] and not r["h"]]
single_h = [r["h"] for r in rows if r["h"] and not r["a"]]
def q(vals, p):
    vals = sorted(v for v in vals if v is not None)
    return vals[int(len(vals) * p)] if vals else None
for name, s in (("alpha-only", single_a), ("HMG-only", single_h)):
    ids = [x["identity"] for x in s]; bs = [x.get("bitscore") for x in s]
    print(f"{name}: n={len(s)} identity p5/p25/p50 = {q(ids,.05):.1f}/{q(ids,.25):.1f}/{q(ids,.5):.1f}"
          f"  bitscore p5/p25/p50 = {q(bs,.05)}/{q(bs,.25)}/{q(bs,.5)}  null bitscore {sum(b is None for b in bs)}")
ta = q([x["identity"] for x in single_a], .05); th = q([x["identity"] for x in single_h], .05)
print("strong threshold = 5th percentile identity of single-gene loci: alpha %.1f, HMG %.1f" % (ta, th))
both = [r for r in rows if r["a"] and r["h"]]
print("both-core loci", len(both))
cat = collections.Counter(); cat_cls = collections.defaultdict(collections.Counter); agree = collections.Counter()
examples = collections.defaultdict(list)
for r in both:
    sa = r["a"]["identity"] >= ta; sh = r["h"]["identity"] >= th
    k = "both strong" if sa and sh else "alpha strong, HMG weak" if sa else "HMG strong, alpha weak" if sh else "both weak"
    cat[k] += 1
    cat_cls[k][r["locus_class"]] += 1
    lab = r["idiomorph"]
    stronger = "MAT1-1" if (r["a"].get("bitscore") or 0) > (r["h"].get("bitscore") or 0) else "MAT1-2"
    agree[(k, lab == stronger)] += 1
    if len(examples[k]) < 4: examples[k].append((r["genome"], round(r["a"]["identity"], 1), round(r["h"]["identity"], 1), lab))
for k, v in cat.most_common():
    ag = agree[(k, True)]; dis = agree[(k, False)]
    print(f"{k:24s} {v:5d}  label = higher-bitscore gene: {ag}, label != higher-bitscore gene: {dis}   loci classes {dict(cat_cls[k])}")
    print("    e.g.", examples[k])
# per-genome: both-strong loci in heterothallic-looking genome? list counts per class for 'both strong' labelled homothallic_candidate
hc = sum(1 for r in both if r["locus_class"] == "homothallic_candidate")
print("both-core loci flagged homothallic_candidate:", hc)
# methods of the weak partner
wm = collections.Counter()
for r in both:
    sa = r["a"]["identity"] >= ta; sh = r["h"]["identity"] >= th
    if sa != sh:
        w = r["h"] if sa else r["a"]
        wm[(w["gene"], w["method"], w["status"])] += 1
print("weak partner (gene, method, status):", wm.most_common(8))

print("\n--- second pass: modelled vs unmodelled, and bitscore gap ---")
def modelled(g): return g["status"] in ("polished_agree", "polished_single", "polished_disagree")
cats = collections.Counter(); gaps = collections.defaultdict(list); lab_dis = collections.Counter(); hc = collections.Counter()
for r in both:
    ma, mh = modelled(r["a"]), modelled(r["h"])
    k = "both modelled" if ma and mh else "alpha modelled only" if ma else "HMG modelled only" if mh else "neither modelled"
    cats[k] += 1
    gap = (r["a"].get("bitscore") or 0) - (r["h"].get("bitscore") or 0)
    gaps[k].append(abs(gap))
    stronger = "MAT1-1" if gap > 0 else "MAT1-2"
    lab_dis[(k, r["idiomorph"] == stronger)] += 1
    if r["locus_class"] == "homothallic_candidate": hc[k] += 1
for k, v in cats.most_common():
    g = sorted(gaps[k])
    print(f"{k:22s} {v:5d}  |bitscore gap| p10/p50/p90 = {g[len(g)//10]:.0f}/{g[len(g)//2]:.0f}/{g[len(g)*9//10]:.0f}  "
          f"gap<20 bits: {sum(x < 20 for x in g)}  label==higher-bitscore {lab_dis[(k, True)]}, differs {lab_dis[(k, False)]}  homothallic_candidate {hc[k]}")
# by class for both modelled
byc = collections.Counter(); 
import re
meta = {}
for row in csv.DictReader(open("/Users/jstajich/projects/MATPredict/results/2026-10-03_ascomycota_v060/loci.tsv"), delimiter="\t"):
    meta[(row["genome"], row["contig"], int(row["start"]))] = row["class_"]
for r in both:
    if modelled(r["a"]) and modelled(r["h"]):
        byc[meta.get((r["genome"], r["contig"], int(r["start"])), "?")] += 1
print("both modelled by class:", byc.most_common(8))
