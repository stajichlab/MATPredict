import csv, sys, statistics, collections
R = "/Users/jstajich/projects/MATPredict/.claude/worktrees/agent-a8e00f1f9146e378f/results/2026-10-18_basidio_scope_only/redPR_skew/"
rows = list(csv.DictReader(open(R + "calls_baseline.tsv"), delimiter="\t"))
called = [r for r in rows if r["label"]]
print("rows", len(rows), "loci", len(called), "genomes with locus", len({r["genome"] for r in called}), "uncalled genomes", sum(1 for r in rows if not r["label"]))
c = collections.Counter(r["label"] for r in called); print("loci", dict(c))
gen = collections.Counter(r["genome"] for r in called if r["label"])
print("\nby genus: genomes(called) A1-loci A2-loci")
for g in sorted({r["genus"] for r in rows}):
    a = collections.Counter(r["label"] for r in called if r["genus"] == g)
    print(f"  {g:20s} genomes={len({r['genome'] for r in rows if r['genus']==g})} A1={a['A1']} A2={a['A2']}")
# per species
sp = collections.defaultdict(lambda: collections.Counter()); spg = collections.defaultdict(set); sptax = collections.defaultdict(set)
for r in rows:
    spg[r["species"]].add(r["genome"]); sptax[r["species"]].add(r["taxid"])
for r in called: sp[r["species"]][r["label"]] += 1
only1 = [s for s, a in sp.items() if a["A1"] and not a["A2"]]; only2 = [s for s, a in sp.items() if a["A2"] and not a["A1"]]; both = [s for s, a in sp.items() if a["A1"] and a["A2"]]
print(f"\nspecies with calls: {len(sp)}; A1 only {len(only1)}, A2 only {len(only2)}, both {len(both)}")
print("species-level (n genomes, A1 loci, A2 loci), species with >=4 genomes or any A1:")
for s, a in sorted(sp.items(), key=lambda t: -len(spg[t[0]])):
    if len(spg[s]) >= 4 or a["A1"]:
        print(f"  {s:42s} genomes={len(spg[s]):3d} taxids={len(sptax[s])} A1={a['A1']:3d} A2={a['A2']:3d}")
# per-taxid dominance:
tx = collections.defaultdict(collections.Counter)
for r in called: tx[r["taxid"]][r["label"]] += 1
print("\ntaxids with calls", len(tx), " taxid A1-only", sum(1 for a in tx.values() if a['A1'] and not a['A2']), " A2-only", sum(1 for a in tx.values() if a['A2'] and not a['A1']), " both", sum(1 for a in tx.values() if a['A1'] and a['A2']))
# margins
def q(v):
    v = sorted(v); n = len(v)
    return f"n={n} min={v[0]:.0f} p10={v[int(n*.1)]:.0f} median={statistics.median(v):.0f} p90={v[int(n*.9)]:.0f} max={v[-1]:.0f}"
for lab in ("A1", "A2"):
    s = [r for r in called if r["label"] == lab]
    win = [float(r["bits_a1" if lab == "A1" else "bits_a2"] or 0) for r in s]
    los = [float(r["bits_a2" if lab == "A1" else "bits_a1"] or 0) for r in s]
    wid = [float(r["id_a1" if lab == "A1" else "id_a2"] or 0) for r in s]
    lid = [float(r["id_a2" if lab == "A1" else "id_a1"] or 0) for r in s]
    mar = [float(r["bits_margin"] or 0) for r in s]
    print(f"\n{lab}: winner bits {q(win)}\n    loser bits {q(los)}\n    winner identity {q(wid)}\n    loser identity {q(lid)}\n    bit margin {q(mar)}")
allm = [float(r["bits_margin"] or 0) for r in called]
for t in (20, 50, 100, 200): print(f"calls with bit margin < {t}: {sum(m < t for m in allm)} of {len(allm)}")
low_id = [r for r in called if float(r["id_a1" if r["label"] == "A1" else "id_a2"] or 0) < 60]
print("calls whose winner identity < 60%:", len(low_id), collections.Counter(r["label"] for r in low_id))
# nearest reference
print("\nnearest (winning) reference record by label:")
ref = collections.Counter((r["label"], r["ref_a1" if r["label"] == "A1" else "ref_a2"]) for r in called)
for k, v in sorted(ref.items()): print("  ", k, v)
# does call genus follow reference genus
def rg(x):
    return {"217517": "Sporobolomyces", "5003": "Sporobolomyces", "5005": "Sporobolomyces"}.get(x.split("_")[0], "Rhodotorula")
same = sum(1 for r in called if rg(r["ref_a1" if r["label"] == "A1" else "ref_a2"]) == r["genus"]); print("winning-ref genus == genome genus:", same, "of", len(called))
# who wins by identity: A1 winners' loser coverage
print("\nloci per genome >1:", [(g, [r['label'] for r in called if r['genome'] == g]) for g, n in gen.items() if n > 1])
