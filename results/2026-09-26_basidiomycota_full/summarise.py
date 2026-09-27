"""Summaries from genomes.tsv / loci.tsv (plain Python, no pandas).

Writes by_order.tsv, misses_by_order.tsv, size_bins.tsv and prints headline
tables. Run with /usr/bin/python3.12 after analyze_full.py.
"""
import csv, collections, statistics, os

HERE = os.path.dirname(os.path.abspath(__file__))
G = list(csv.DictReader(open(os.path.join(HERE, "genomes.tsv")), delimiter="\t"))
L = list(csv.DictReader(open(os.path.join(HERE, "loci.tsv")), delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def med(v):
    v = [x for x in v if x is not None]
    return round(statistics.median(v)) if v else ""


def miss_cause(g):
    """Why an uncalled genome is uncalled. Order of tests matters."""
    if g["status"] == "no_report":
        return "no_report (timeout/fail)"
    n50, contigs = num(g["n50"]), num(g["contigs"])
    sup_pol = num(g["sup_max_polished"])
    sup_genes = num(g["sup_max_genes"])
    if g["n_suppressed"] not in ("", "0") and sup_pol is not None and sup_pol >= 1:
        cause = "bar: withheld locus with >=1 modelled gene"
    elif g["n_suppressed"] not in ("", "0"):
        cause = "bar: withheld locus, 0 modelled"
    else:
        cause = "nothing localised (reference gap)"
    if n50 is not None and n50 < 20000:
        cause += " + poor assembly (N50<20kb)"
    return cause


# ---- by subphylum / order
by_sub = collections.defaultdict(lambda: collections.Counter())
by_ord = collections.defaultdict(lambda: collections.Counter())
walls = collections.defaultdict(list)
for g in G:
    for key, d in ((g["subphylum"] or "?", by_sub), ((g["subphylum"] or "?", g["order"] or "?"), by_ord)):
        d[key]["n"] += 1
        d[key][g["status"]] += 1
        d[key]["route_" + (g["routing"] or "none")] += 1
    walls[(g["subphylum"] or "?", g["order"] or "?")].append(num(g["wall_s"]))

print("== by subphylum")
for k, c in sorted(by_sub.items(), key=lambda kv: -kv[1]["n"]):
    print(f"{k:28s} n={c['n']:5d} called={c['called']:5d} ({100*c['called']/c['n']:.1f}%) uncalled={c['uncalled']} no_report={c['no_report']}")

rows = []
for k, c in sorted(by_ord.items(), key=lambda kv: -kv[1]["n"]):
    w = walls[k]
    rows.append(dict(subphylum=k[0], order=k[1], n=c["n"], called=c["called"],
                     pct=round(100 * c["called"] / c["n"], 1), uncalled=c["uncalled"], no_report=c["no_report"],
                     lineage=c["route_lineage"], phylum_fallback=c["route_phylum_fallback"],
                     not_searched=c["route_not_searched"], wall_median_s=med(w),
                     wall_max_s=max([x for x in w if x is not None], default="")))
with open(os.path.join(HERE, "by_order.tsv"), "w", newline="") as fo:
    wr = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); wr.writeheader(); wr.writerows(rows)
print("\n== by order (n >= 10)")
for r in rows:
    if r["n"] >= 10:
        print(f"{r['subphylum'][:18]:18s} {r['order'][:22]:22s} n={r['n']:4d} called={r['called']:4d} {r['pct']:5.1f}%  lin={r['lineage']} fb={r['phylum_fallback']} ns={r['not_searched']} wall med/max={r['wall_median_s']}/{r['wall_max_s']}")

# ---- families / confidence / class / idiomorph
print("\n== loci by family x confidence")
fc = collections.Counter((l["family_called"], l["confidence"]) for l in L)
for f in sorted({l["family_called"] for l in L}):
    print(f"  {f:28s}", {c: fc[(f, c)] for c in ("high", "medium", "low") if fc[(f, c)]},
          dict(collections.Counter(l["locus_class"] for l in L if l["family_called"] == f)),
          "idiomorph:", dict(collections.Counter(l["idiomorph"] for l in L if l["family_called"] == f).most_common(6)))
print("detection_pass:", dict(collections.Counter(l["detection_pass"] for l in L)))
fam_by_sub = collections.Counter((l["subphylum"], l["family_called"]) for l in L)
print("loci by subphylum x family:", dict(fam_by_sub))
per_genome = collections.Counter(l["genome"] for l in L)
print("loci per called genome:", dict(sorted(collections.Counter(per_genome.values()).items())))
fams_per_genome = collections.Counter(g["families_called"] for g in G if g["status"] == "called")
print("family combinations:", dict(fams_per_genome.most_common(8)))

# ---- flags
tot = lambda k: sum(int(num(g[k]) or 0) for g in G)
print("\nflags: homothallic", tot("homothallic"), "unverified", tot("unverified"),
      "genomes with assembly_gap_at_locus", sum(1 for g in G if (num(g["gap_at_locus"]) or 0) > 0),
      "idiomorph_unmodelled loci", tot("idiomorph_unmodelled"), "suppressed_flank_carried", tot("suppressed_flank_carried"),
      "zygosity values", dict(collections.Counter(g["zygosity"] for g in G)))

# ---- misses
mc = collections.defaultdict(collections.Counter)
for g in G:
    if g["status"] != "called":
        mc[(g["subphylum"] or "?", g["order"] or "?")][miss_cause(g)] += 1
causes = sorted({c for v in mc.values() for c in v})
with open(os.path.join(HERE, "misses_by_order.tsv"), "w", newline="") as fo:
    wr = csv.writer(fo, delimiter="\t"); wr.writerow(["subphylum", "order", "uncalled"] + causes)
    for k, v in sorted(mc.items(), key=lambda kv: -sum(kv[1].values())):
        wr.writerow([k[0], k[1], sum(v.values())] + [v[c] for c in causes])
allc = collections.Counter()
for v in mc.values():
    allc.update(v)
print("\n== miss causes, all:", dict(allc.most_common()))
print("== miss causes, top orders:")
for k, v in sorted(mc.items(), key=lambda kv: -sum(kv[1].values()))[:12]:
    print(f"  {k[1][:24]:24s} {sum(v.values()):4d}", dict(v.most_common()))

# ---- re-run effect proxies
cap = [g for g in G if (num(g["capped_clusters"]) or 0) > 0]
print("\n== genomes with >=1 polish-capped cluster:", len(cap), "of which uncalled:", sum(1 for g in cap if g["status"] != "called"))
relax = [g for g in G if g["status"] != "called" and (num(g["nd_best_fraction"]) or 0) < 0.5
         and (num(g["sup_max_genes"]) or 0) >= 2]
print("uncalled genomes with best fraction < 0.5 and a withheld cluster of >=2 genes (relaxed-pass upper bound):", len(relax))

# ---- size bins vs wall
bins = [(0, 50e6), (50e6, 100e6), (100e6, 200e6), (200e6, 500e6), (500e6, 1e12)]
names = ["<50 Mb", "50-100 Mb", "100-200 Mb", "200-500 Mb", ">500 Mb"]
out = []
for (lo, hi), nm in zip(bins, names):
    gg = [g for g in G if num(g["size_bp"]) is not None and lo <= num(g["size_bp"]) < hi]
    w = sorted(x for x in (num(g["wall_s"]) for g in gg) if x is not None)
    to = sum(1 for g in gg if g["status"] == "no_report" or (num(g["wall_s"]) or 0) >= 3590)
    p90 = w[int(0.9 * (len(w) - 1))] if w else ""
    out.append(dict(bin=nm, genomes=len(gg), wall_median_s=med(w), wall_p90_s=p90, wall_max_s=w[-1] if w else "",
                    timeouts_or_no_report=to,
                    orders="; ".join(f"{o}:{n}" for o, n in collections.Counter(g["order"] for g in gg).most_common(4))))
with open(os.path.join(HERE, "size_bins.tsv"), "w", newline="") as fo:
    wr = csv.DictWriter(fo, fieldnames=list(out[0]), delimiter="\t"); wr.writeheader(); wr.writerows(out)
print("\n== size bins")
for r in out:
    print(r)
nosize = [g["genome"] for g in G if num(g["size_bp"]) is None]
print("genomes without asm_stats size:", len(nosize))
EOF_MARK = None
