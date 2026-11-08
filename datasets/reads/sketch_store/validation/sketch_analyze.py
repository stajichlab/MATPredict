"""Can filtered read sketches (scaled=100, k21/k31, abund>=2) of ~30x subsamples support PCoA and MAT/target-gene calls in A. fumigatus?
Inputs in W: sigs/<strain>.{A,B}.min2.sig.zip, strains.tsv, mat/<strain>.tsv, snp_gt.tsv.gz (bcftools query GT matrix, header = samples),
targets/*.fasta (optional; one record per target) and targets/truth_overlap.tsv (optional). Writes W/analysis/."""
import csv, glob, gzip, math, os, sys, collections, itertools
import numpy as np
import sourmash
from sourmash import MinHash

W = sys.argv[1]
OUT = os.path.join(W, "analysis_" + os.environ.get("MODE", "solid")); os.makedirs(OUT, exist_ok=True)
strains = [l.split("\t")[0] for l in open(os.path.join(W, "strains.tsv")) if l.strip()]
kind = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in open(os.path.join(W, "strains.tsv")) if l.strip()}
KS = (21, 31)
MAXH = 2 ** 64


def load(path, k):
    if os.path.getsize(path) < 1000: return None
    try:
        for ss in sourmash.load_file_as_signatures(path, ksize=k):
            return dict(ss.minhash.hashes) or None
    except Exception:
        return None
    return None


sig = {k: {} for k in KS}
for s in strains:
    for h in "AB":
        p = os.path.join(W, "sigs", f"{s}.{h}.min2.sig.zip")
        if os.path.exists(p):
            for k in KS:
                d = load(p, k)
                if d: sig[k][(s, h)] = d
have = [s for s in strains if (s, "A") in sig[21]]
print(f"strains with sketches: {len(have)} of {len(strains)}")
rep = open(os.path.join(OUT, "REPORT.txt"), "w")
def say(*a):
    print(*a); print(*a, file=rep)


class HS:
    """hash set of one sample: all hashes (abund>=2 already) and the solid ones (abund >= max(3, SOLID x median abund))"""
    def __init__(self, d, scaled):
        lim = MAXH // scaled
        d = {h: n for h, n in d.items() if h < lim}
        self.all = set(d)
        med = float(np.median(list(d.values()))) if d else 0
        thr = max(3, SOLID * med)
        self.solid = {h for h, n in d.items() if n >= thr}
    def __len__(self): return len(self.all)

SOLID = float(os.environ.get("SOLID", "0.25"))
MODE = os.environ.get("MODE", "solid")


def ani(a, b, k):
    """containment ANI, mean of both directions (C^(1/k)). MODE=plain: |A∩B|/|A| over all hashes.
    MODE=solid: |solid(A) ∩ all(B)| / |solid(A)|, so error k-mers (low abundance) leave the denominator and coverage gaps in B
    (k-mers seen once or twice) still count as shared."""
    if not a.all or not b.all: return float("nan")
    if MODE == "plain":
        i = len(a.all & b.all); return 0.5 * ((i / len(a.all)) ** (1 / k) + (i / len(b.all)) ** (1 / k))
    ca = len(a.solid & b.all) / max(1, len(a.solid)); cb = len(b.solid & a.all) / max(1, len(b.solid))
    return 0.5 * (ca ** (1 / k) + cb ** (1 / k))


# ---- 1. SNP distances
snp = None
gtp = os.path.join(W, "snp_gt.tsv.gz")
if os.path.exists(gtp):
    with gzip.open(gtp, "rt") as fh:
        names = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            b = line.rstrip("\n").encode()
            v = np.frombuffer(b[::2], dtype=np.uint8).astype(np.int8) - 48 if len(b) == 2 * len(names) - 1 else \
                np.array([(-1 if g[0] == "." else int(g[0])) for g in line.rstrip("\n").split("\t")], dtype=np.int8)
            v[v < 0] = -1
            rows.append(v)
    G = np.vstack(rows)
    idx = {n: i for i, n in enumerate(names)}
    X1 = (G == 1).astype(np.float32); X0 = (G == 0).astype(np.float32)
    DIFF = X1.T @ X0 + X0.T @ X1
    snp = {}
    for a, b in itertools.combinations([s for s in have if s in idx], 2):
        snp[(a, b)] = int(DIFF[idx[a], idx[b]])
    say(f"SNP matrix: {G.shape[0]} sites x {G.shape[1]} samples; pairs with SNP distance: {len(snp)}")

# ---- 2. distances, ceiling, agreement with SNPs, per k and scaled
summary = []
best = None
for k in KS:
    for scaled in (100, 500, 1000):
        A = {s: HS(sig[k][(s, "A")], scaled) for s in have}
        B = {s: HS(sig[k][(s, "B")], scaled) for s in have if (s, "B") in sig[k]}
        ceil = {s: ani(A[s], B[s], k) for s in B}
        pairs = {(a, b): ani(A[a], A[b], k) for a, b in itertools.combinations(have, 2)}
        cv = sorted(v for v in ceil.values() if not math.isnan(v))
        pv = sorted(pairs.values())
        lowc = cv[0] if cv else float("nan"); medc = cv[len(cv) // 2] if cv else float("nan")
        reach = sum(v >= lowc for v in pv)
        rho = float("nan")
        if snp:
            common = [p for p in pairs if p in snp]
            x = np.array([1 - pairs[p] for p in common]); y = np.array([snp[p] for p in common])
            rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y))
            rho = float(np.corrcoef(rx, ry)[0, 1])
        nh = int(np.median([len(A[s]) for s in have]))
        summary.append((k, scaled, nh, lowc, medc, pv[len(pv) // 2], pv[int(0.99 * len(pv))], reach, len(pv), rho))
        if k == 21 and scaled == 100:
            best = (A, B, ceil, pairs)
        with open(os.path.join(OUT, f"pairs_k{k}_s{scaled}.tsv"), "w") as fh:
            fh.write("a\tb\tani\tsnps\n")
            for (a, b), v in pairs.items():
                fh.write(f"{a}\t{b}\t{v:.6f}\t{snp.get((a, b), '') if snp else ''}\n")
        with open(os.path.join(OUT, f"ceiling_k{k}_s{scaled}.tsv"), "w") as fh:
            fh.write("strain\tA_vs_B_ani\thashes_A\thashes_B\n")
            for s, v in ceil.items():
                fh.write(f"{s}\t{v:.6f}\t{len(A[s])}\t{len(B[s])}\n")
say("\nk  scaled  median_hashes  sameStrain_min  sameStrain_median  between_median  between_p99  between>=sameStrain_min  spearman(sketch dist, SNPs)")
for r in summary:
    say(f"{r[0]}  {r[1]:5}  {r[2]:8}  {r[3]:.5f}  {r[4]:.5f}  {r[5]:.5f}  {r[6]:.5f}  {r[7]}/{r[8]}  {r[9]:.3f}")

# lab-strain families: known same-strain groups
fam = {"Af293": ["AF293", "Af293", "AF293-VKBR1", "AF293-VKBR3", "AF293-VKWH2"],
       "CEA10": ["CEA10", "CEA10-2", "CEA10-MCH", "CEA10_CramerSM", "A1163", "akuBKU80"],
       "A1159": ["A1159", "A1159-T0-2"], "ATCC46645": ["ATCC_46645", "ATCC_46645.DMC"]}
A, B, ceil, pairs = best
say("\nKnown same-strain groups (k21 scaled100): within-group sketch ANI and SNPs; rank of the closest non-member")
for f, mem in fam.items():
    mem = [m for m in mem if m in A]
    for a, b in itertools.combinations(mem, 2):
        p = (a, b) if (a, b) in pairs else (b, a)
        say(f"  {f}: {a} vs {b}: ani {pairs[p]:.5f}  SNPs {snp.get(p, snp.get((p[1], p[0]), '')) if snp else ''}")
    if mem and len(mem) < len(A):
        outside = [v for (a, b), v in pairs.items() if (a in mem) != (b in mem) and (a in mem or b in mem)]
        say(f"  {f}: closest strain outside the group: ani {max(outside):.5f}")

# ---- 3. PCoA of sketch distances vs PCA of SNPs
def pcoa(D, n=3):
    m = D.shape[0]; J = np.eye(m) - 1 / m
    Bm = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(Bm); o = np.argsort(w)[::-1]; w, v = w[o], v[:, o]
    return v[:, :n] * np.sqrt(np.maximum(w[:n], 0)), w[:n] / w[w > 0].sum()
pop = [s for s in have if kind[s] != "lab" or s in ("AF293", "CEA10")]
ix = {s: i for i, s in enumerate(pop)}
D = np.zeros((len(pop), len(pop)))
for (a, b), v in pairs.items():
    if a in ix and b in ix: D[ix[a], ix[b]] = D[ix[b], ix[a]] = 1 - v
X, ve = pcoa(D)
say(f"\nPCoA of sketch distances ({len(pop)} strains): variance explained by axes 1-3: {', '.join(f'{x:.3f}' for x in ve)}")
if snp:
    Ds = np.zeros_like(D)
    for (a, b), v in snp.items():
        if a in ix and b in ix: Ds[ix[a], ix[b]] = Ds[ix[b], ix[a]] = v
    Y, ves = pcoa(Ds)
    say(f"PCoA of SNP distances: variance explained 1-3: {', '.join(f'{x:.3f}' for x in ves)}")
    for i in range(3):
        say(f"  axis {i+1}: |r| sketch vs SNP = {abs(np.corrcoef(X[:, i], Y[:, i])[0, 1]):.3f}")
    # Procrustes-like: correlation of the two 2-D configurations after optimal rotation
    from numpy.linalg import svd
    def norm(Z): Z = Z - Z.mean(0); return Z / np.linalg.norm(Z)
    a2, b2 = norm(X[:, :2]), norm(Y[:, :2]); u, s_, vt = svd(a2.T @ b2)
    say(f"  Procrustes correlation (axes 1-2): {s_.sum():.3f}")
    with open(os.path.join(OUT, "pcoa.tsv"), "w") as fh:
        fh.write("strain\tsk1\tsk2\tsk3\tsnp1\tsnp2\tsnp3\n")
        for s in pop:
            i = ix[s]; fh.write(f"{s}\t" + "\t".join(f"{z:.5f}" for z in list(X[i]) + list(Y[i])) + "\n")
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
        grp = {l.split("\t")[0]: l.split("\t")[1].split("/")[0] for l in open(os.path.join(W, "strains.tsv"))}
        gs = sorted({grp[s] for s in pop}); cmap = plt.get_cmap("tab20")
        for j, (Z, t) in enumerate(((X, "read sketches (k21, scaled=100)"), (Y, "VCF SNP distances"))):
            for gi, g in enumerate(gs):
                pts = [ix[s] for s in pop if grp[s] == g]
                ax[j].scatter(Z[pts, 0], Z[pts, 1], s=14, color=cmap(gi % 20), label=g)
            ax[j].set_title(t); ax[j].set_xlabel("PCo1"); ax[j].set_ylabel("PCo2")
        ax[1].legend(fontsize=6, ncol=2, loc="best"); fig.tight_layout(); fig.savefig(os.path.join(OUT, "pcoa.png"), dpi=130)
    except Exception as e:
        say("plot failed:", e)

# ---- 4. target screens: MAT idiomorph cores and any targets/*.fasta
def read_fasta(p):
    recs, name = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"): name = line[1:].split()[0]; recs[name] = []
        elif name: recs[name].append(line)
    return {n: "".join(v) for n, v in recs.items()}

def hashes_of(seq, k, scaled=100):
    mh = MinHash(n=0, ksize=k, scaled=scaled); mh.add_sequence(seq.upper(), force=True); return set(mh.hashes)

panel = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave1/plan/aspergillus_fumigatus/panel/L1"
r1, r2 = read_fasta(panel + "/aspergillus_MAT1-1.fasta"), read_fasta(panel + "/aspergillus_MAT1-2.fasta")
mat_calls = {}
for s in have:
    p = os.path.join(W, "mat", s + ".tsv")
    if os.path.exists(p):
        rows_ = list(csv.DictReader(open(p), delimiter="\t"))
        if rows_: mat_calls[s] = rows_[0]
for k in KS:
    c1 = set().union(*[hashes_of(v, k) for n, v in r1.items() if "core" in n])
    c2 = set().union(*[hashes_of(v, k) for n, v in r2.items() if "core" in n])
    fl = set().union(*[hashes_of(v, k) for n, v in {**r1, **r2}.items() if "flank" in n])
    u1, u2 = c1 - c2 - fl, c2 - c1 - fl
    say(f"\nMAT, k={k}: idiomorph-specific hashes at scaled=100: MAT1-1 {len(u1)}, MAT1-2 {len(u2)}; flank hashes {len(fl)}")
    tab = collections.Counter(); rows_out = []
    for s in have:
        hs = sig[k][(s, "A")]
        n1, n2 = len(u1 & hs.keys()), len(u2 & hs.keys())
        f1, f2 = n1 / max(1, len(u1)), n2 / max(1, len(u2))
        call = "MAT1-1" if f1 >= 0.5 and f2 < 0.1 else "MAT1-2" if f2 >= 0.5 and f1 < 0.1 else "both" if f1 >= 0.5 and f2 >= 0.5 else "none/unclear"
        rc = mat_calls.get(s, {}).get("call", "?")
        tab[(rc, call)] += 1
        rows_out.append((s, n1, len(u1), n2, len(u2), call, rc))
    with open(os.path.join(OUT, f"mat_sketch_k{k}.tsv"), "w") as fh:
        fh.write("strain\thits_MAT1-1\tof\thits_MAT1-2\tof\tsketch_call\treads_type_call\n")
        for r in rows_out: fh.write("\t".join(map(str, r)) + "\n")
    say("  reads-type call -> sketch call: " + ", ".join(f"{a}->{b}: {n}" for (a, b), n in sorted(tab.items(), key=lambda x: -x[1])))
    say("  hits per strain (MAT1-1 / MAT1-2), first 8: " + "; ".join(f"{r[0]} {r[1]}/{r[3]}" for r in rows_out[:8]))

tf = sorted(glob.glob(os.path.join(W, "targets", "*.fa*")))
if tf:
    truth = {}
    tp = os.path.join(W, "targets", "truth_overlap.tsv")
    if os.path.exists(tp):
        for r in csv.DictReader(open(tp), delimiter="\t"):
            truth[r.get("strain_global_name") or list(r.values())[0]] = r
    k = 21
    tg = {}
    for p in tf:
        tg.update(read_fasta(p))
    say(f"\nTarget screen, k={k} scaled=100: {len(tg)} targets from {', '.join(os.path.basename(p) for p in tf)}")
    with open(os.path.join(OUT, "targets_k21.tsv"), "w") as fh:
        fh.write("strain\ttarget\ttarget_len\ttarget_hashes\thits\tfraction\n")
        per = collections.defaultdict(list)
        for t, seq in tg.items():
            th = hashes_of(seq, k)
            for s in have:
                hs = sig[k][(s, "A")]
                n = len(th & hs.keys())
                fh.write(f"{s}\t{t}\t{len(seq)}\t{len(th)}\t{n}\t{n / max(1, len(th)):.3f}\n")
                per[t].append(n / max(1, len(th)))
        for t, v in per.items():
            v = sorted(v)
            say(f"  {t[:50]:50} len {len(tg[t]):7}  hashes {len(hashes_of(tg[t], k)):6}  strains with fraction>=0.8: {sum(x >= 0.8 for x in v)}, 0.2-0.8: {sum(0.2 <= x < 0.8 for x in v)}, <0.2: {sum(x < 0.2 for x in v)}")
    # agreement with truth: fraction of target hashes found vs truth presence; best threshold per target and overall at 0.5
    frac = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(OUT, "targets_k21.tsv")), delimiter="\t"):
        frac[r["target"]][r["strain"]] = float(r["fraction"])
    say("\nTarget vs truth (strains with truth and a sketch); call present if fraction >= 0.5. Columns: TP FN FP TN; truth-positive fractions; truth-negative max")
    for t in tg:
        col = "bafA" if t.lower().startswith("bafa") or "afu5g14915" in t.lower() else t
        tp = fn = fp = tn = 0; posf, negf = [], []
        for s, f in frac[t].items():
            tv = truth.get(s, {}).get(col)
            if tv not in ("0", "1"): continue
            if tv == "1": posf.append(f); tp += f >= 0.5; fn += f < 0.5
            else: negf.append(f); fp += f >= 0.5; tn += f < 0.5
        if posf or negf:
            say(f"  {t[:24]:24} TP {tp:3} FN {fn:3} FP {fp:3} TN {tn:3} | positives: {' '.join(f'{x:.2f}' for x in sorted(posf)[:12])} | negatives max {max(negf) if negf else float('nan'):.2f}, median {np.median(negf) if negf else float('nan'):.2f}")
rep.close()
