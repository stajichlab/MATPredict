"""Floors on the strongest core hit of each flank-carried locus.
Current: E <= 1e-5 genome-wide. Candidates: bitscore >= 30..50; option (b): E rescaled to a fixed
search space of 1e8 residues (E_fixed = E * 1e8 / (2 * genome_nt), six frames of L/3)."""
import csv, collections
rows = list(csv.DictReader(open("scored.tsv"), delimiter="\t"))
def f(x, d=0.0):
    try: return float(x)
    except: return d
def grp(r):
    s = r["source"]
    return {"audit_cap6": "Ascomycota (cap6 audit)", "serinales_882aa01": "Serinales MTL", "mucoro_4457c3a": "Mucoromycota",
            "early_div_634dda4": "Mucoro-group 634dda4", "basidio_full_ad1f865": "Basidiomycota HD",
            "dothideo_recheck_9932c7c": "Dothideomycetes"}.get(s, s)
rules = {"E<=1e-5 (current)": lambda r: r["evalue"] != "" and f(r["evalue"], 99) <= 1e-5}
for b in (30, 35, 40, 45, 50): rules[f"bits>={b}"] = (lambda b: lambda r: f(r["bits"]) >= b)(b)
def efix(r):
    L = f(r["genome_len"]); e = f(r["evalue"], 99)
    return e * 1e8 / (2 * L) if L and r["evalue"] != "" else 99
rules["E_fixed1e8<=1e-5"] = lambda r: efix(r) <= 1e-5
out = open("tradeoffs.txt", "w")
def p(*a): print(*a); print(*a, file=out)
real = [r for r in rows if r["label"] == "real"]; noise = [r for r in rows if r["label"] == "noise"]
best = [r for r in rows if r["label"] == "unknown" and r["rank_genomewide"] == "1"]
notbest = [r for r in rows if r["label"] == "unknown" and r["rank_genomewide"] not in ("1",)]
p(f"loci {len(rows)}; real {len(real)}; independent noise {len(noise)}; unlabelled best-locus {len(best)}; unlabelled not-best {len(notbest)}")
p("\nreal loci:"); [p("  ", r["source"], r["species"], r["best_gene"], r["bits"], r["evalue"], r["genome_len"]) for r in real]
p("\nrule                real_kept  noise_kept  bestlocus_kept  notbest_kept  total_pass")
for k, fn in rules.items():
    p(f"{k:20s} {sum(map(fn,real)):3d}/{len(real):<5d} {sum(map(fn,noise)):3d}/{len(noise):<6d} {sum(map(fn,best)):4d}/{len(best):<8d} {sum(map(fn,notbest)):4d}/{len(notbest):<7d} {sum(map(fn,rows)):4d}/{len(rows)}")
p("\npass counts per group")
G = sorted({grp(r) for r in rows})
p("rule                " + "  ".join(f"{g[:22]:>22s}" for g in G))
for k, fn in rules.items():
    p(f"{k:20s}" + "  ".join(f"{sum(fn(r) for r in rows if grp(r)==g):>10d}/{sum(1 for r in rows if grp(r)==g):<11d}" for g in G))
p("\nbitscore distribution by label (min / median / max)")
import statistics as st
for lab, xs in (("real", real), ("noise", noise), ("best-locus", best), ("not-best", notbest)):
    v = sorted(f(r["bits"]) for r in xs)
    if v: p(f"  {lab:10s} n={len(v):3d}  {v[0]:.1f} / {st.median(v):.1f} / {v[-1]:.1f}")
p("\nnoise loci passing bits>=40:")
for r in noise:
    if f(r["bits"]) >= 40: p("  ", r["species"], r["genome"], r["best_gene"], r["bits"], r["evalue"], r["pid"])
