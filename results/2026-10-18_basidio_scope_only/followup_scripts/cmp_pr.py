import csv, sys, collections
R = "/Users/jstajich/projects/MATPredict/.claude/worktrees/agent-a8e00f1f9146e378f/results/2026-10-18_basidio_scope_only/redPR_skew/"
def load(f):
    d = collections.defaultdict(list)
    for r in csv.DictReader(open(R + f), delimiter="\t"):
        if r["label"]: d[r["genome"]].append(r["label"])
        else: d.setdefault(r["genome"], [])
    return d
base = load("calls_baseline.tsv")
for name, f in (("E1 one A1+one A2 per genus (4 refs)", "calls_E1_one_per_allele_per_genus.tsv"), ("E3 A2 scarce", "calls_E3_a2_scarce.tsv"), ("E4 A1 scarce", "calls_E4_a1_scarce.tsv")):
    try: v = load(f)
    except FileNotFoundError: continue
    gs = sorted(set(v) & set(base))
    tot = collections.Counter(l for g in gs for l in v[g]); tb = collections.Counter(l for g in gs for l in base[g])
    same = sum(sorted(v[g]) == sorted(base[g]) for g in gs)
    changed = [(g, base[g], v[g]) for g in gs if sorted(v[g]) != sorted(base[g])]
    lost = sum(1 for g in gs if base[g] and not v[g]); gained = sum(1 for g in gs if not base[g] and v[g])
    flip = sum(1 for g in gs if base[g] and v[g] and sorted(base[g]) != sorted(v[g]) and len(base[g]) == len(v[g]) == 1)
    print(f"\n{name}: genomes {len(gs)}; baseline loci {dict(tb)}; variant loci {dict(tot)}; identical genome calls {same}; lost {lost}; gained {gained}; A1<->A2 flips {flip}")
    for c in changed[:12]: print("   ", c)
