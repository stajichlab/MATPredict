#!/usr/bin/env python3
"""Queries over a sketch store (STORE/sigs/<sample_id>.sig.zip: k=21,k=31 scaled=100 FracMinHash with abundances, k-mers seen
once dropped; built by sketch_one.sh). Each subcommand reads the sketches and writes STORE/query/<name>/.

  qc        per sample: hashes, solid hashes, median abundance (a coverage proxy) and all/solid ratio (high = mixture or
            contamination: the test set's normal range was 1.0-1.1, mixtures 1.4-1.8)
  distances all-against-all "solid containment" ANI (k=21, downsampled to scaled=1000: same answers as scaled=100 on the test set),
            outliers, PCoA, clone groups (average linkage) and a FastME tree
  screen    presence of target sequences (FASTA; one record per target): fraction of each target's specific hashes found and their
            median abundance relative to the sample's median (a copy-number estimate); optional truth table
  mat       MAT idiomorph calls from idiomorph-specific hashes (MATPredict reads panel), beside the reads-type calls in STORE/mat

Solid containment: a hash counts in sample a's denominator only if its abundance in a is at least max(3, 0.25 x a's median), and as
shared if b has it at all, so sequencing errors leave the denominator and coverage gaps in b do not count as differences.
ANI(a,b) = mean over both directions of containment^(1/k). Validated 2026-11-07 on 136 strains against SNP calls
(results/2026-11-07_sketch_pcoa_mat): same-strain median 0.99996, between-strain p99 0.999; PCoA Procrustes vs SNPs 0.63-0.70.
"""
import argparse, collections, csv, glob, gzip, json, math, os, subprocess, sys
from multiprocessing import Pool
import numpy as np

MAXH = 2 ** 64
SOLID = 0.25


def sample_files(store):
    fs = sorted(glob.glob(os.path.join(store, "sigs", "*.sig.zip")))
    return [(os.path.basename(f)[:-len(".sig.zip")], f) for f in fs]


def load(path, k=21):
    import sourmash
    for ss in sourmash.load_file_as_signatures(path, ksize=k):
        return dict(ss.minhash.hashes)
    return {}


def solid_thr(d):
    med = float(np.median(list(d.values()))) if d else 0.0
    return med, max(3.0, SOLID * med)


# ---------------------------------------------------------------- qc
def _qc(arg):
    sid, path = arg
    d = load(path)
    med, thr = solid_thr(d)
    ns = sum(1 for v in d.values() if v >= thr)
    return sid, len(d), ns, med


def cmd_qc(a):
    out = os.path.join(a.store, "query", "qc"); os.makedirs(out, exist_ok=True)
    with Pool(a.threads) as p:
        rows = p.map(_qc, sample_files(a.store), chunksize=8)
    stats = {}
    for f in glob.glob(os.path.join(a.store, "stats", "*.tsv")):
        for line in open(f):
            x = line.rstrip("\n").split("\t"); stats[x[0]] = x
    ratios = [r[1] / max(1, r[2]) for r in rows]
    med_ratio = float(np.median(ratios))
    with open(os.path.join(out, "qc.tsv"), "w") as fh:
        fh.write("sample_id\thashes_k21_s100\tsolid_hashes\tmedian_abund\tall_over_solid\tbases_used\tcoverage_used\tflag\n")
        for (sid, n, ns, med), r in zip(rows, ratios):
            st = stats.get(sid, [""] * 5)
            flag = ";".join(f for f, c in (("mixture?", r > 1.3 * med_ratio), ("low_coverage", med < 8), ("few_hashes", ns < 0.7 * np.median([x[2] for x in rows]))) if c)
            fh.write(f"{sid}\t{n}\t{ns}\t{med:.1f}\t{r:.3f}\t{st[3]}\t{st[4]}\t{flag}\n")
    print(f"qc: {len(rows)} samples; median all/solid {med_ratio:.3f}; written {out}/qc.tsv")


# ---------------------------------------------------------------- distances
def _down(arg):
    sid, path, scaled = arg
    d = load(path)
    lim = MAXH // scaled
    d = {h: v for h, v in d.items() if h < lim}
    med, thr = solid_thr(d)
    allh = np.fromiter(d.keys(), dtype=np.uint64, count=len(d))
    sol = np.fromiter((h for h, v in d.items() if v >= thr), dtype=np.uint64)
    return sid, allh, sol


def cmd_distances(a):
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    out = os.path.join(a.store, "query", "distances"); os.makedirs(out, exist_ok=True)
    files = sample_files(a.store)
    with Pool(a.threads) as p:
        res = p.map(_down, [(s, f, a.scaled) for s, f in files], chunksize=8)
    ids = [r[0] for r in res]; n = len(ids); k = 21
    # columns: hashes seen in at least two samples (a hash in one sample adds nothing to any intersection)
    cnt = collections.Counter()
    for _, al, _ in res: cnt.update(al.tolist())
    cols = {h: i for i, h in enumerate(h for h, c in cnt.items() if c >= 2)}
    H = len(cols)
    print(f"distances: {n} samples, {H} shared hashes at scaled={a.scaled}", flush=True)
    A = np.zeros((n, H), dtype=np.float32); S = np.zeros((n, H), dtype=np.float32)
    nsol = np.zeros(n)
    for i, (_, al, so) in enumerate(res):
        A[i, [cols[h] for h in al.tolist() if h in cols]] = 1
        S[i, [cols[h] for h in so.tolist() if h in cols]] = 1
        nsol[i] = len(so)
    I = S @ A.T                                   # I[a,b] = |solid(a) ∩ all(b)|
    del A, S
    C = I / np.maximum(nsol, 1)[:, None]
    np.clip(C, 1e-12, 1, out=C)
    C = C ** (1 / k)
    ANI = 0.5 * (C + C.T); np.fill_diagonal(ANI, 1.0)
    np.save(os.path.join(out, "ani_k21.npy"), ANI.astype(np.float32))
    with open(os.path.join(out, "samples.txt"), "w") as fh: fh.write("\n".join(ids) + "\n")
    # outliers: mean ANI to all others far below the rest
    mean = (ANI.sum(1) - 1) / (n - 1); m = np.median(mean); mad = np.median(np.abs(mean - m))
    outl = mean < m - 8 * max(mad, 1e-6)
    keep = np.where(~outl)[0]
    print(f"outliers (mean ANI < median - 8 MAD): {int(outl.sum())}", flush=True)
    D = 1 - ANI[np.ix_(keep, keep)].astype(np.float64)
    # PCoA (classical MDS), top 10 axes
    from scipy.linalg import eigh
    m_ = len(keep); J = np.eye(m_) - 1.0 / m_
    B = -0.5 * J @ (D ** 2) @ J
    w, v = eigh(B, subset_by_index=[m_ - 10, m_ - 1]); o = np.argsort(w)[::-1]; w, v = w[o], v[:, o]
    X = v * np.sqrt(np.maximum(w, 0)); tot = np.trace(B)
    with open(os.path.join(out, "pcoa.tsv"), "w") as fh:
        fh.write("sample_id\t" + "\t".join(f"PCo{i+1}" for i in range(10)) + "\n")
        for j, i in enumerate(keep): fh.write(ids[i] + "\t" + "\t".join(f"{x:.6g}" for x in X[j]) + "\n")
    # clone groups: average linkage on the kept samples at several identity thresholds
    Z = linkage(squareform(D, checks=False), method="average")
    thr = [0.9999, 0.9995, 0.999, 0.998]
    labs = {t: fcluster(Z, 1 - t, criterion="distance") for t in thr}
    with open(os.path.join(out, "samples_qc.tsv"), "w") as fh:
        fh.write("sample_id\tmean_ani\toutlier\t" + "\t".join(f"clone_{t}" for t in thr) + "\n")
        pos = {i: j for j, i in enumerate(keep)}
        for i, s in enumerate(ids):
            cl = "\t".join(str(labs[t][pos[i]]) if i in pos else "" for t in thr)
            fh.write(f"{s}\t{mean[i]:.6f}\t{int(outl[i])}\t{cl}\n")
    # tree: FastME (balanced minimum evolution, NNI+SPR) on 1-ANI, names as indices (mapping in tree_names.tsv)
    ph = os.path.join(out, "dist.phy")
    with open(ph, "w") as fh:
        fh.write(f"{m_}\n")
        for j in range(m_): fh.write(f"S{j:05d} " + " ".join(f"{x:.7f}" for x in D[j]) + "\n")
    with open(os.path.join(out, "tree_names.tsv"), "w") as fh:
        for j, i in enumerate(keep): fh.write(f"S{j:05d}\t{ids[i]}\n")
    fastme = a.fastme or "/opt/linux/rocky/8.x/x86_64/pkgs/fastme/2.1.6.3/bin/fastme"
    tree = os.path.join(out, "tree_fastme.nwk")
    r = subprocess.run([fastme, "-i", ph, "-o", tree, "-m", "B", "-n", "-s", "-T", str(a.threads)], capture_output=True, text=True)
    if r.returncode == 0 and os.path.exists(tree):
        nm = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(out, "tree_names.tsv")))
        t = open(tree).read()
        import re
        t = re.sub(r"S\d{5}", lambda mm: nm.get(mm.group(0), mm.group(0)), t)
        open(os.path.join(out, "tree_fastme.named.nwk"), "w").write(t)
    else:
        print("fastme failed:", r.stderr[-500:])
    iu = np.triu_indices(m_, 1)
    print(f"PCoA: variance share of axes 1-3 {', '.join(f'{x / tot:.3f}' for x in w[:3])}; between-sample ANI median "
          f"{np.median(1 - D[iu]):.5f}, p99 {np.quantile(1 - D[iu], .99):.5f}; clone groups at 0.9995: {len(set(labs[0.9995]))}")


# ---------------------------------------------------------------- target hits (shared by screen and mat)
def read_fasta(p):
    recs, name = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"): name = line[1:].split()[0]; recs[name] = []
        elif name: recs[name].append(line)
    return {n: "".join(v) for n, v in recs.items()}


def hashes_of(seq, k=21, scaled=100):
    from sourmash import MinHash
    mh = MinHash(n=0, ksize=k, scaled=scaled); mh.add_sequence(seq.upper(), force=True); return set(mh.hashes)


_PROBE = None
def _hits(arg):
    sid, path = arg
    d = load(path)
    med = float(np.median(list(d.values()))) if d else 0.0
    return sid, med, {h: d[h] for h in _PROBE if h in d}


def probe_all(store, probe, threads):
    global _PROBE
    _PROBE = probe
    with Pool(threads) as p:
        return p.map(_hits, sample_files(store), chunksize=8)


def solo_hashes(spec, hit_sets, max_solo=0.01, min_samples=10):
    """Hashes of one target's specific set that are not specific in practice: seen ALONE (none of the target's other specific
    hashes in the sample) in more than max(min_samples, max_solo x samples) samples. Such a hash is shared with sequence outside
    the target (another allele, a paralog, a repeat): in the A. fumigatus store one of bafA's three hashes was alone in 778 isolates
    that lack the gene and its Starship. Returns {hash: samples where it is alone}; a target with one specific hash is left as is."""
    if len(spec) < 2 or not max_solo:
        return {}
    solo = collections.Counter()
    for hits in hit_sets:
        got = [h for h in spec if h in hits]
        if len(got) == 1:
            solo[got[0]] += 1
    lim = max(min_samples, max_solo * len(hit_sets))
    drop = {h: c for h, c in solo.items() if c > lim}
    return drop if len(drop) < len(spec) else {}      # never drop every hash of a target


def cmd_screen(a):
    out = os.path.join(a.store, "query", "screen_" + a.name); os.makedirs(out, exist_ok=True)
    tg = {}
    for f in a.targets: tg.update(read_fasta(f))
    H = {t: hashes_of(s) for t, s in tg.items()}
    res = probe_all(a.store, set().union(*H.values()), a.threads)
    n = len(res)
    present_in = collections.Counter(h for _, _, hits in res for h in hits)
    core = {h for h, c in present_in.items() if c >= a.core_fraction * n}
    linked = [tuple(x.split(",")) for x in (a.linked or [])]
    spec = {}
    for t in tg:
        others = set().union(*[H[u] for u in H if u != t and (t, u) not in linked and (u, t) not in linked])
        spec[t] = H[t] - others - core
    dropped = {}
    for t in tg:
        d = solo_hashes(spec[t], [hits for _, _, hits in res], a.max_solo)
        if d:
            spec[t] = spec[t] - set(d); dropped[t] = d
    with open(os.path.join(out, "dropped_hashes.tsv"), "w") as fh:
        fh.write("target\thash\tsamples_alone\n")
        for t, d in dropped.items():
            for h, c in sorted(d.items(), key=lambda x: -x[1]):
                fh.write(f"{t}\t{h}\t{c}\n")
    truth = {}
    if a.truth:
        for r in csv.DictReader(open(a.truth), delimiter="\t"):
            truth[r[a.truth_id]] = r
    rows = []
    with open(os.path.join(out, "screen.tsv"), "w") as fh:
        fh.write("sample_id\ttarget\ttarget_len\ttarget_hashes\tspecific_hashes\thits\tfraction\trel_abund\tcall\n")
        for sid, med, hits in res:
            for t in tg:
                sp = spec[t]; got = [hits[h] for h in sp if h in hits]
                f = len(got) / len(sp) if sp else float("nan")
                rel = (float(np.median(got)) / med) if got and med else 0.0
                call = "present" if f >= a.min_fraction else "partial" if f >= 0.2 else "absent"
                rows.append((sid, t, f, call))
                fh.write(f"{sid}\t{t}\t{len(tg[t])}\t{len(H[t])}\t{len(sp)}\t{len(got)}\t{f:.3f}\t{rel:.2f}\t{call}\n")
    with open(os.path.join(out, "summary.txt"), "w") as fh:
        for t in tg:
            c = collections.Counter(call for s, tt, f, call in rows if tt == t)
            line = f"{t[:40]:40} len {len(tg[t]):7}  specific hashes {len(spec[t]):6}  present {c['present']}  partial {c['partial']}  absent {c['absent']}"
            if truth:
                col = a.truth_map.get(t, t) if a.truth_map else t
                tp = fn = fp = tn = 0
                for s, tt, f, call in rows:
                    if tt != t: continue
                    v = truth.get(s, {}).get(col)
                    if v not in ("0", "1"): continue
                    p = call == "present"
                    tp += p and v == "1"; fn += (not p) and v == "1"; fp += p and v == "0"; tn += (not p) and v == "0"
                line += f"  | vs truth TP {tp} FN {fn} FP {fp} TN {tn}"
            if t in dropped:
                line += f"  | dropped {len(dropped[t])} hash(es) seen alone in {max(dropped[t].values())} samples"
            print(line); fh.write(line + "\n")


def cmd_mat(a):
    out = os.path.join(a.store, "query", "mat"); os.makedirs(out, exist_ok=True)
    r1 = read_fasta(os.path.join(a.panel, f"{a.prefix}_MAT1-1.fasta")); r2 = read_fasta(os.path.join(a.panel, f"{a.prefix}_MAT1-2.fasta"))
    c1 = set().union(*[hashes_of(v) for n, v in r1.items() if "core" in n]); c2 = set().union(*[hashes_of(v) for n, v in r2.items() if "core" in n])
    fl = set().union(*[hashes_of(v) for n, v in {**r1, **r2}.items() if "flank" in n])
    u1, u2 = c1 - c2 - fl, c2 - c1 - fl
    print(f"mat: idiomorph-specific hashes MAT1-1 {len(u1)}, MAT1-2 {len(u2)}")
    res = probe_all(a.store, u1 | u2, a.threads)
    man = {r["sample_id"]: r for r in csv.DictReader(open(os.path.join(a.store, "manifest.tsv")), delimiter="\t")}
    tab = collections.Counter()
    with open(os.path.join(out, "mat_calls.tsv"), "w") as fh:
        fh.write("sample_id\tsource\tstudy_accession\tcountry\tsketch_hits_MAT1-1\tof\tsketch_hits_MAT1-2\tof\tsketch_call\treads_type_call\treads_type_flags\tscreen_status\n")
        for sid, med, hits in res:
            n1 = sum(1 for h in u1 if h in hits); n2 = sum(1 for h in u2 if h in hits)
            f1, f2 = n1 / max(1, len(u1)), n2 / max(1, len(u2))
            call = ("both" if f1 >= .5 and f2 >= .5 else "MAT1-1" if f1 >= .5 and f2 < .2 else "MAT1-2" if f2 >= .5 and f1 < .2 else "none/unclear")
            m = man.get(sid, {}); key = m.get("run_accession") or sid
            rt, flags, scr = "", "", ""
            p = os.path.join(a.store, "mat", f"{key}.{a.prefix}.tsv")
            if os.path.exists(p):
                rr = list(csv.DictReader(open(p), delimiter="\t"))
                if rr: rt, flags = rr[-1].get("call", ""), rr[-1].get("flags", "")
            sp = os.path.join(a.store, "mat", f"{key}.screen.json")
            if os.path.exists(sp):
                try: scr = json.load(open(sp)).get("status", "")
                except Exception: scr = "?"
            tab[(rt or "-", call)] += 1
            fh.write(f"{sid}\t{m.get('source', '')}\t{m.get('study_accession', '')}\t{m.get('country', '')}\t{n1}\t{len(u1)}\t{n2}\t{len(u2)}\t{call}\t{rt}\t{flags}\t{scr}\n")
    s = "reads-type -> sketch: " + ", ".join(f"{x}->{y}: {c}" for (x, y), c in sorted(tab.items(), key=lambda z: -z[1]))
    print(s); open(os.path.join(out, "summary.txt"), "w").write(s + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--store", required=True); ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("qc")
    d = sub.add_parser("distances"); d.add_argument("--scaled", type=int, default=1000); d.add_argument("--fastme")
    s = sub.add_parser("screen"); s.add_argument("--name", required=True); s.add_argument("--targets", nargs="+", required=True)
    s.add_argument("--min-fraction", type=float, default=0.5); s.add_argument("--core-fraction", type=float, default=0.95)
    s.add_argument("--linked", nargs="*", help="TARGET_A,TARGET_B pairs whose shared hashes are kept (e.g. a gene and the element that carries it)")
    s.add_argument("--max-solo", type=float, default=0.01, help="drop a specific hash seen alone (none of its target's other specific hashes "
                   "present) in more than this fraction of samples (and more than 10); 0 keeps every hash. Dropped hashes: dropped_hashes.tsv")
    s.add_argument("--truth"); s.add_argument("--truth-id", default="sample_id")
    s.add_argument("--truth-map", type=lambda x: dict(kv.split("=", 1) for kv in x.split(";")), help="target=truth_column;...")
    m = sub.add_parser("mat"); m.add_argument("--panel", required=True); m.add_argument("--prefix", required=True)
    a = ap.parse_args()
    {"qc": cmd_qc, "distances": cmd_distances, "screen": cmd_screen, "mat": cmd_mat}[a.cmd](a)


if __name__ == "__main__":
    main()
