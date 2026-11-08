"""Follow-up: outliers, PCoA without them, and target scoring on target-specific hashes."""
import csv, os, sys, collections, itertools
import numpy as np
sys.argv = [sys.argv[0], "."]
exec(open("sketch_analyze.py").read().split("# ---- 1. SNP distances")[0])   # loads sig, have, kind, HS, ani, load, MODE=solid
W = "."
OUT = "analysis_followup"; os.makedirs(OUT, exist_ok=True)
rep = open(os.path.join(OUT, "REPORT.txt"), "w")
def say(*a):
    print(*a); print(*a, file=rep)
k = 21
A = {s: HS(sig[k][(s, "A")], 100) for s in have}
B = {s: HS(sig[k][(s, "B")], 100) for s in have if (s, "B") in sig[k]}
pairs = {}
for r in csv.DictReader(open("analysis_solid/pairs_k21_s100.tsv"), delimiter="\t"):
    pairs[(r["a"], r["b"])] = (float(r["ani"]), int(r["snps"]) if r["snps"] else None)
def d(a, b): return pairs.get((a, b), pairs.get((b, a)))
ceil = {s: ani(A[s], B[s], k) for s in B}
grp = {l.split("\t")[0]: l.split("\t")[1].split("/")[0] for l in open("strains.tsv")}
meanani = {s: np.mean([d(s, t)[0] for t in have if t != s]) for s in have}
stats = {l.split("\t")[0] + l.split("\t")[1]: l.rstrip().split("\t") for f in os.listdir("stats") for l in open("stats/" + f)}
say("Lowest split-sample ceilings and lowest mean ANI to others (k21 scaled100 solid):")
for s in sorted(ceil, key=ceil.get)[:8]:
    st = stats.get(s + "A", [""] * 7)
    say(f"  ceiling {ceil[s]:.5f}  mean ani {meanani[s]:.5f}  {s} ({grp[s]})  hashes A {len(A[s])} solid {len(A[s].solid)}  readlen {st[2]}")
say("  ...")
for s in sorted(meanani, key=meanani.get)[:8]:
    say(f"  mean ani {meanani[s]:.5f}  ceiling {ceil.get(s, float('nan')):.5f}  {s} ({grp[s]})  hashes {len(A[s])}")
cv = sorted(v for v in ceil.values())
say(f"ceiling quantiles: min {cv[0]:.5f} p05 {cv[int(.05*len(cv))]:.5f} median {cv[len(cv)//2]:.5f}")

# outliers: mean ANI more than 5 MAD below the median
m = np.median(list(meanani.values())); mad = np.median([abs(v - m) for v in meanani.values()])
out = {s for s, v in meanani.items() if v < m - 8 * max(mad, 1e-6)}
say(f"\noutliers (mean ANI < median - 8 MAD): {sorted(out)}")
keep = [s for s in have if s not in out and (kind[s] != "lab" or s in ("AF293", "CEA10"))]
def pcoa(D, n=3):
    mm = D.shape[0]; J = np.eye(mm) - 1 / mm
    Bm = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(Bm); o = np.argsort(w)[::-1]; w, v = w[o], v[:, o]
    return v[:, :n] * np.sqrt(np.maximum(w[:n], 0)), w[:n] / w[w > 0].sum()
ix = {s: i for i, s in enumerate(keep)}
Dk = np.zeros((len(keep),) * 2); Ds = np.zeros_like(Dk)
for a, b in itertools.combinations(keep, 2):
    v, sn = d(a, b); Dk[ix[a], ix[b]] = Dk[ix[b], ix[a]] = 1 - v; Ds[ix[a], ix[b]] = Ds[ix[b], ix[a]] = sn
X, ve = pcoa(Dk); Y, ves = pcoa(Ds)
iu = np.triu_indices(len(keep), 1)
rk = lambda z: np.argsort(np.argsort(z))
say(f"\nwithout outliers and lab copies: {len(keep)} strains, {len(iu[0])} pairs")
say(f"  spearman(sketch distance, SNPs) = {np.corrcoef(rk(Dk[iu]), rk(Ds[iu]))[0,1]:.3f};  pearson = {np.corrcoef(Dk[iu], Ds[iu])[0,1]:.3f}")
say(f"  sketch PCoA variance 1-3: {', '.join(f'{x:.3f}' for x in ve)};  SNP PCoA: {', '.join(f'{x:.3f}' for x in ves)}")
for i in range(3): say(f"  axis {i+1}: |r| = {abs(np.corrcoef(X[:, i], Y[:, i])[0,1]):.3f}")
def norm(Z): Z = Z - Z.mean(0); return Z / np.linalg.norm(Z)
for n in (2, 3):
    u, s_, vt = np.linalg.svd(norm(X[:, :n]).T @ norm(Y[:, :n])); say(f"  Procrustes correlation (axes 1-{n}): {s_.sum():.3f}")
# nearest neighbour agreement: is each strain's closest strain by sketch among its 3 closest by SNPs?
nn_ok = 0
for s in keep:
    i = ix[s]; dk = Dk[i].copy(); dk[i] = np.inf; ds = Ds[i].copy(); ds[i] = np.inf
    nn_ok += int(np.argmin(dk) in set(np.argsort(ds)[:3]))
say(f"  nearest neighbour by sketch is among the 3 nearest by SNPs: {nn_ok}/{len(keep)}")
with open(os.path.join(OUT, "pcoa.tsv"), "w") as fh:
    fh.write("strain\tgroup\tsk1\tsk2\tsk3\tsnp1\tsnp2\tsnp3\n")
    for s in keep: fh.write(f"{s}\t{grp[s]}\t" + "\t".join(f"{z:.6f}" for z in list(X[ix[s]]) + list(Y[ix[s]])) + "\n")
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.8)); gs = sorted({grp[s] for s in keep}); cmap = plt.get_cmap("tab20")
    for j, (Z, t) in enumerate(((X, "read sketches (k21, scaled=100, solid containment)"), (Y, "SNP distances (VCF)"))):
        for gi, g in enumerate(gs):
            pts = [ix[s] for s in keep if grp[s] == g]; ax[j].scatter(Z[pts, 0], Z[pts, 1], s=14, color=cmap(gi % 20), label=g)
        ax[j].set_title(t, fontsize=9); ax[j].set_xlabel("PCo1"); ax[j].set_ylabel("PCo2")
    ax[1].legend(fontsize=6, ncol=2); fig.tight_layout(); fig.savefig(os.path.join(OUT, "pcoa.png"), dpi=130)
except Exception as e:
    say("plot failed", e)

# ---- targets on specific hashes
def read_fasta(p):
    recs, name = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"): name = line[1:].split()[0]; recs[name] = []
        elif name: recs[name].append(line)
    return {n: "".join(v) for n, v in recs.items()}
from sourmash import MinHash
def hashes_of(seq, k, scaled=100):
    mh = MinHash(n=0, ksize=k, scaled=scaled); mh.add_sequence(seq.upper(), force=True); return set(mh.hashes)
tg = {**read_fasta("targets/bafA.fa"), **read_fasta("targets/starships.fa")}
H = {t: hashes_of(s, k) for t, s in tg.items()}
# hashes in (nearly) every strain are core genome, not element-specific
cnt = collections.Counter(h for s in have for h in A[s].all if any(h in H[t] for t in H))
core = {h for h, n in cnt.items() if n >= 0.95 * len(have)}
truth = {r["strain_global_name"]: r for r in csv.DictReader(open("targets/truth_overlap.tsv"), delimiter="\t")}
say(f"\nTargets on SPECIFIC hashes (minus every other target's hashes and minus hashes in >=95% of strains); present if fraction >= 0.5")
say("target                    all  specific | TP FN FP TN | positives (lowest 5) | negatives max")
rows = []
for t in tg:
    others = set().union(*[H[u] for u in H if u != t and not (t.startswith("bafA") and u == "Nebuchadnezzar-h1") and not (t == "Nebuchadnezzar-h1" and u.startswith("bafA"))])
    sp = H[t] - others - core
    col = "bafA" if t.startswith("bafA") else t
    tp = fn = fp = tn = 0; pf, nf = [], []
    for s in have:
        f = len(sp & A[s].all) / len(sp) if sp else float("nan")
        rows.append((s, t, len(H[t]), len(sp), f))
        tv = truth.get(s, {}).get(col)
        if tv == "1": pf.append(f); tp += f >= .5; fn += f < .5
        elif tv == "0": nf.append(f); fp += f >= .5; tn += f < .5
    say(f"{t[:24]:24} {len(H[t]):5} {len(sp):8} | {tp:3} {fn:3} {fp:3} {tn:3} | {' '.join(f'{x:.2f}' for x in sorted(pf)[:5])} | {max(nf) if nf else float('nan'):.2f}")
with open(os.path.join(OUT, "targets_specific_k21.tsv"), "w") as fh:
    fh.write("strain\ttarget\ttarget_hashes\tspecific_hashes\tfraction\n")
    for r in rows: fh.write("\t".join(map(str, r[:4])) + f"\t{r[4]:.3f}\n")
# disagreements for bafA / Nebuchadnezzar-h1
say("\nbafA / Nebuchadnezzar-h1 disagreements with truth (strain, truth bafA, source, sketch bafA frac, sketch Neb-h1 frac):")
fr = {(r[0], r[1]): r[4] for r in rows}
bt = [t for t in tg if t.startswith("bafA")][0]
for s in have:
    tv = truth.get(s, {}).get("bafA")
    if tv in ("0", "1") and (fr[(s, bt)] >= .5) != (tv == "1"):
        say(f"  {s}: truth {tv} ({truth[s].get('bafA_source')}), sketch bafA {fr[(s, bt)]:.2f}, Neb-h1 {fr[(s, 'Nebuchadnezzar-h1')]:.2f}")
rep.close()
