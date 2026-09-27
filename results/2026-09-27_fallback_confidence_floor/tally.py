"""Count phylum_fallback (and, separately, explicit_phylum) calls below each floor."""
import csv, collections
rows = list(csv.DictReader(open("calls.tsv"), delimiter="\t"))
FLOORS = (30, 35, 40)
out = open("tally.txt", "w")


def p(*a):
    print(*a); print(*a, file=out)


p("routing modes by run (calls):")
for run in sorted({r["run"] for r in rows}):
    c = collections.Counter(r["routing"] for r in rows if r["run"] == run)
    p(f"  {run:20s} {dict(c)}")

for routing in ("phylum_fallback", "explicit_phylum"):
    p(f"\n=== {routing} ===")
    sel = [r for r in rows if r["routing"] == routing]
    for run in sorted({r["run"] for r in sel}):
        R = [r for r in sel if r["run"] == run]
        for conf in ("high", "medium"):
            C = [r for r in R if r["confidence"] == conf]
            line = f"  {run:20s} {conf:6s} n={len(C):4d}"
            for var in ("best_core_any", "best_core_modelled"):
                cells = []
                for fl in FLOORS:
                    k = sum(1 for r in C if r[var] != "" and float(r[var]) < fl)
                    cells.append(f"<{fl}:{k}")
                nomod = sum(1 for r in C if r[var] == "")
                line += f" | {var.split('_')[-1]} {' '.join(cells)} none:{nomod}"
            p(line)

p("\n=== phylum_fallback HIGH calls that would drop to medium, by run/family/order (best_core_any) ===")
sel = [r for r in rows if r["routing"] == "phylum_fallback" and r["confidence"] == "high"]
g = collections.defaultdict(lambda: [0, 0, 0, 0])
for r in sel:
    key = (r["run"], r["family"], r["order"])
    g[key][3] += 1
    if r["best_core_any"] != "":
        v = float(r["best_core_any"])
        for i, fl in enumerate(FLOORS):
            if v < fl:
                g[key][i] += 1
p(f"  {'run':18s} {'family':24s} {'order':22s} <30 <35 <40 high_total")
for k, v in sorted(g.items(), key=lambda kv: -kv[1][2]):
    if v[2]:
        p(f"  {k[0]:18s} {k[1]:24s} {k[2]:22s} {v[0]:3d} {v[1]:3d} {v[2]:3d} {v[3]:4d}")
