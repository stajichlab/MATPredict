"""One arm: determinism, full-build HMM stats, gate threshold on 189 negatives,
genus LOO on 108 positives (training + Zygo) as full/hmgbox/50-90 aa windows.
Usage: run_arm.py ARM [TRIM]   TRIM in {none, clipkit, gap50, hmgbox}"""
import json, random, subprocess, sys, tempfile, time
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import *  # noqa

arm = sys.argv[1]; trimname = sys.argv[2] if len(sys.argv) > 2 else "none"
RES = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/results")
F3 = RES / "2026-09-28_validation_f3_f4"
AA = pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile(str(RES / "2026-09-26_sexMP_phylogeny/PF00505.hmm")) as fh:
    PF = fh.read()
CLIPKIT = "/opt/linux/rocky/8.x/x86_64/pkgs/clipkit/1.3.0/bin/clipkit"


def env_of(seqs):
    block = pyhmmer.easel.DigitalSequenceBlock(AA, [pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(AA) for k, v in seqs.items()])
    out = {}
    for hit in pyhmmer.plan7.Pipeline(AA, E=10, domE=10).search_hmm(PF, block):
        n = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
        d = max(hit.domains, key=lambda d: d.score)
        out[n] = (d.env_from - 1, d.env_to)
    return out


def write_cols(aln, keep, out):
    a = cb.read_fasta(aln)
    cb.write_fasta(out, {k: "".join(v[i] for i in keep) for k, v in a.items()})
    return out


def trim(aln: Path):
    if trimname == "clipkit":
        out = aln.with_suffix(".ck.afa")
        subprocess.run([CLIPKIT, str(aln), "-m", "kpic-gappy", "-o", str(out)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return out
    a = cb.read_fasta(aln); L = len(next(iter(a.values()))); n = len(a)
    if trimname == "gap50":
        keep = [i for i in range(L) if sum(v[i] == "-" for v in a.values()) <= 0.5 * n]
        return write_cols(aln, keep, aln.with_suffix(".g50.afa"))
    if trimname == "hmgbox":
        env = env_of({k: v.replace("-", "") for k, v in a.items()})
        lo, hi = L, 0
        for k, (s, e) in env.items():
            pos = [i for i, c in enumerate(a[k]) if c != "-"]
            lo = min(lo, pos[s]); hi = max(hi, pos[min(e, len(pos)) - 1])
        keep = list(range(max(0, lo - 10), min(L, hi + 11)))
        return write_cols(aln, keep, aln.with_suffix(".hmg.afa"))
    return aln


patch(arm, None if trimname == "none" else trim)
t0 = time.time()
fam = family(); CL = CODE / "db/Mucoromycota/classifiers/MAT"
rows = cb.training_set(CODE / "db", fam, CL)
res = dict(arm=arm, trim=trimname)

# 1 determinism: three untrimmed full alignments per gene, and three HMM files
det = {}
for gene in ("sexP", "sexM"):
    seqs = {r["id"]: r["sequence"] for r in rows if r["gene"] == gene}
    hs = set(); hh = set()
    for i in range(3):
        w = Path(tempfile.mkdtemp())
        fa = w / "in.faa"; cb.write_fasta(fa, cb._sorted(seqs))
        hs.add(sha(align(arm, fa, w / "a.afa")))
        hmm, _ = cb.build_hmm(seqs, gene, w)
        buf = w / "h.hmm"; hmm.write(open(buf, "wb")); hh.add(sha(buf))
    det[gene] = dict(alignments_identical=len(hs) == 1, hmms_identical=len(hh) == 1)
res["determinism"] = det

# 2 full build stats + gate
w = Path(tempfile.mkdtemp()); full = {}; stats = {}
bg = pyhmmer.plan7.Background(AA)
for gene in ("sexP", "sexM"):
    seqs = {r["id"]: r["sequence"] for r in rows if r["gene"] == gene}
    hmm, aln = cb.build_hmm(seqs, gene, w)
    full[gene] = hmm
    ncol = len(next(iter(cb.read_fasta(aln).values())))
    stats[gene] = dict(match_states=hmm.M, aln_columns=ncol, mean_match_relative_entropy=round(hmm.mean_match_relative_entropy(bg), 3))
res["hmm_stats"] = stats
neg = cb.read_fasta(CL / "paralog_negatives.faa")
negsc = {g: cb.score(full[g], neg) for g in full}
best = [max(negsc[g][k] for g in full) for k in neg]
res["gate_threshold"] = cb.gate_threshold(best)

# 3 genus LOO with fragments (same fragment generator/seed as F3)
zy = cb.read_fasta(RES / "2026-09-26_sexMP_hmm/zygo_locus_proteins.faa")
zdom = {l.split()[0] for l in open(RES / "2026-09-26_sexMP_hmm/zygo_pf.domtbl") if not l.startswith("#")}
zyrows = [dict(id=k, gene="sexP" if k.split("|")[2] == "Plus" else "sexM", genus=k.split("|")[1].split("_")[0], source="zygo", sequence=v) for k, v in zy.items() if k in zdom]
rng = random.Random(11)
allpos = rows + zyrows
env = env_of({r["id"]: r["sequence"] for r in allpos})
frags = {}
for r in sorted(allpos, key=lambda r: r["id"]):
    k, s = r["id"], r["sequence"]
    frags[(k, "full")] = s
    if k in env:
        a, b = env[k]
        frags[(k, "hmgbox")] = s[max(0, a - 5):min(len(s), b + 5)]
        mid = (a + b) // 2
        for i in range(3):
            L = rng.randint(50, 90)
            lo = max(0, min(mid - rng.randint(5, L - 5), len(s) - L))
            if len(s) >= 50:
                frags[(k, f"win{i}")] = s[lo:lo + L]
out = []
for g in sorted({r["genus"] for r in allpos}):
    train = [r for r in rows if r["genus"] != g]
    hm = {}
    for gene in ("sexM", "sexP"):
        seqs = {r["id"]: r["sequence"] for r in train if r["gene"] == gene}
        hm[gene], _ = cb.build_hmm(seqs, f"loo_{gene}", w)
    held = {r["id"]: r for r in allpos if r["genus"] == g}
    flat = {f"{k}||{t}": v for (k, t), v in frags.items() if k in held and len(v) >= 20}
    sc = {gene: cb.score(h, flat) for gene, h in hm.items()}
    for fk in flat:
        k, t = fk.split("||"); r = held[k]
        own = sc[r["gene"]][fk]; oth = sc["sexM" if r["gene"] == "sexP" else "sexP"][fk]
        out.append(dict(id=k, source=r["source"], frag="win" if t.startswith("win") else t, own=own, margin=own - oth))
by = defaultdict(list)
for o in out: by[o["frag"]].append(o)
summ = {}
for f, L in by.items():
    m = [o["margin"] for o in L]
    summ[f] = dict(n=len(L), wrong=sum(x <= 0 for x in m), called_at25=sum(x >= 25 for x in m), worst_correct=round(min(x for x in m if x > 0), 1) if any(x > 0 for x in m) else None,
                   median_margin=round(sorted(m)[len(m) // 2], 1))
fullpos = [o for o in by["full"]]
summ["full_true_at_or_above_gate"] = sum(o["own"] >= res["gate_threshold"] for o in fullpos)
summ["full_training_only"] = dict(n=sum(o["source"] != "zygo" for o in fullpos), wrong=sum(o["margin"] <= 0 for o in fullpos if o["source"] != "zygo"),
                                  worst_correct=round(min(o["margin"] for o in fullpos if o["source"] != "zygo" and o["margin"] > 0), 1))
z = [o for o in fullpos if o["source"] == "zygo"]
summ["zygo_full"] = dict(n=len(z), wrong=sum(o["margin"] <= 0 for o in z), min_margin=round(min(o["margin"] for o in z), 1), median_own=round(sorted(o["own"] for o in z)[len(z) // 2], 1))
res["loo"] = summ
res["runtime_s"] = round(time.time() - t0, 1)
tag = f"{arm}__{trimname}"
json.dump(res, open(R / f"res_{tag}.json", "w"), indent=1)
with open(R / f"loo_{tag}.tsv", "w") as fo:
    fo.write("id\tsource\tfrag\town\tmargin\n")
    for o in out: fo.write(f"{o['id']}\t{o['source']}\t{o['frag']}\t{o['own']:.1f}\t{o['margin']:.1f}\n")
print(json.dumps(res))
