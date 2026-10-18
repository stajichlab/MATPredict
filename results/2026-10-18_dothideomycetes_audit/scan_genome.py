"""Genome-wide, annotation-free search for alpha-box and HMG-box ORFs (six-frame translation, stop-to-stop fragments >= 50 aa)
with Pfam (MATalpha_HMGbox PF04769, HMG-box clan) and the curated alpha-box / HMG-box profiles, independent of the
pipeline's own search and polishing.  Answers: does a genome that is called one idiomorph also carry the other idiomorph's
core gene somewhere (homothallic or retroposed fragment), and is there an alpha-box ORF at all behind a MAT1-1 call.

    scan_genome.py --list genomes.tsv --lib LIB --db DB --pfam-subset pfam_subset.hmm --out scan.tsv
genomes.tsv: ASM_ACCESSION column (BFD id).  Standard genetic code only (fragments are scored as protein, so a
non-standard code only changes stops).  Output: one row per genome with the best hits and their positions."""
import argparse, csv, gzip, sys, tempfile, time
from pathlib import Path
import numpy as np

CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"  # TCAG order
MAP = np.full(256, 4, dtype=np.uint8)
for i, ch in enumerate(b"TCAG"):
    MAP[ch] = i
TABLE = np.array(list(CODE) + ["X"], dtype="<U1")
COMP = bytes.maketrans(b"ACGTacgtN", b"TGCAtgcaN")


def read_all(path):
    out, name, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf).upper()
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper()
    return out


def six_frames(seq):
    """yield (strand, frame, aa_string) with positions recoverable: aa index i -> nt = frame + 3*i (+ strand handling)."""
    for strand in ("+", "-"):
        s = seq if strand == "+" else seq.encode().translate(COMP).decode()[::-1]
        arr = np.frombuffer(s.encode(), dtype=np.uint8)
        codes = MAP[arr]
        for f in range(3):
            n = (len(codes) - f) // 3
            if n <= 0:
                continue
            c = codes[f:f + 3 * n].reshape(n, 3).astype(np.int16)
            bad = (c == 4).any(axis=1)
            idx = c[:, 0] * 16 + c[:, 1] * 4 + c[:, 2]
            idx[bad] = 64
            yield strand, f, "".join(TABLE[idx])


def orfs(seq, minlen=50):
    L = len(seq)
    for strand, f, aa in six_frames(seq):
        pos = 0
        for frag0 in aa.split("*"):
            # a stop-free stretch longer than HMMER's 100,000-residue limit (N runs) is cut into overlapping windows
            pieces = [(0, frag0)] if len(frag0) <= 50000 else [(o, frag0[o:o + 30000]) for o in range(0, len(frag0), 29000)]
            for off, frag in pieces:
                if len(frag) >= minlen:
                    a0, a1 = pos + off, pos + off + len(frag)
                    if strand == "+":
                        s, e = f + 3 * a0 + 1, f + 3 * a1
                    else:
                        s, e = L - (f + 3 * a1) + 1, L - (f + 3 * a0)
                    yield strand, s, e, frag
            pos += len(frag0) + 1


def best_hits(hmms, frags, alpha, T):
    """hmmsearch each HMM against the fragments; return {name: [(score, frag_index)...]} for hits >= T (GA when T is None)."""
    import pyhmmer
    from pyhmmer.easel import TextSequence
    seqs = [TextSequence(name=str(i).encode(), sequence=f[3]).digitize(alpha) for i, f in enumerate(frags)]
    res = {}
    for name, hmm in hmms.items():
        out = []
        kw = {"bit_cutoffs": "gathering"} if T is None else {"T": T, "incE": 1e9, "E": 1e9}
        for hits in pyhmmer.hmmsearch(hmm, seqs, cpus=1, **kw):
            for h in hits:
                if T is not None or h.included:
                    out.append((h.score, int(h.name.decode() if isinstance(h.name, bytes) else h.name)))
        res[name] = sorted(out, reverse=True)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", required=True); ap.add_argument("--lib", required=True)
    ap.add_argument("--db", required=True, type=Path); ap.add_argument("--pfam-subset", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    import pyhmmer
    from pyhmmer.easel import Alphabet
    alpha = Alphabet.amino()
    pf = {}
    for h in pyhmmer.plan7.HMMFile(a.pfam_subset):
        pf[h.name.decode() if isinstance(h.name, bytes) else h.name] = h
    want_pf = {k: pf[k] for k in ("MATalpha_HMGbox", "HMG_box", "HMG_box_2") if k in pf}
    prof = {}
    d = a.db / "Ascomycota" / "classifiers" / "MAT_core"
    for g in ("MAT1-1-1", "MAT1-2-1"):
        with pyhmmer.plan7.HMMFile(str(d / f"{g}.hmm")) as fh:
            prof[g] = fh.read()
    asm = [r["ASM_ACCESSION"] for r in csv.DictReader(open(a.list), delimiter="\t")]
    if a.limit:
        asm = asm[:a.limit]
    cols = ["genome", "n_orfs", "alpha_pfam_best", "alpha_pfam_n", "alpha_pfam_where", "alpha_prof_best", "alpha_prof_n40",
            "alpha_prof_where", "hmg_pfam_n", "hmg_prof_best", "hmg_prof_n40", "hmg_prof_where", "alpha_hits40", "hmg_hits40",
            "hmg_pfam_hits", "secs"]
    t0 = time.time()
    done = set()
    if Path(a.out).exists():  # resume: keep genomes already scanned
        done = {r["genome"] for r in csv.DictReader(open(a.out), delimiter="\t")}
    with open(a.out, "a" if done else "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        if not done:
            w.writeheader()
        for i, g in enumerate(asm):
            if g in done:
                continue
            t1 = time.time()
            try:
                contigs = read_all(Path(a.lib) / f"{g}.fa.gz")
            except Exception as exc:  # noqa: BLE001
                w.writerow(dict(genome=g, n_orfs=-1)); continue
            frags = []
            for cname, s in contigs.items():
                for strand, st, en, aa in orfs(s):
                    frags.append((cname, st, en, aa, strand))
            frags2 = [(c, s_, e_, aa) for c, s_, e_, aa, _ in frags]
            # fragments pre-filter by the (fast) alpha/HMG Pfam models; then profile-score only what passes a low bar
            hp = best_hits(want_pf, frags2, alpha, None)
            hq = best_hits(prof, frags2, alpha, 15.0)
            def where(lst):
                if not lst:
                    return ""
                f = frags[lst[0][1]]
                return f"{f[0]}:{f[1]}-{f[2]}{f[4]}"
            ap_, aq, hq_ = hp.get("MATalpha_HMGbox", []), hq["MAT1-1-1"], hq["MAT1-2-1"]
            hp_n = len({i for k in ("HMG_box", "HMG_box_2") for _, i in hp.get(k, [])})
            def hits40(lst, T=40.0):
                out = []
                for sc_, i in lst:
                    if sc_ >= T and len(out) < 12:
                        f = frags[i]
                        out.append(f"{f[0]}:{f[1]}-{f[2]}{f[4]}:{sc_:.0f}:{len(f[3])}")
                return ";".join(out)
            hpf = sorted({(round(sc_, 1), i) for k in ("HMG_box", "HMG_box_2") for sc_, i in hp.get(k, [])}, reverse=True)
            w.writerow(dict(genome=g, alpha_hits40=hits40(aq), hmg_hits40=hits40(hq_), hmg_pfam_hits=hits40(hpf, 0.0),
                            n_orfs=len(frags), alpha_pfam_best=round(ap_[0][0], 1) if ap_ else "",
                            alpha_pfam_n=len(ap_), alpha_pfam_where=where(ap_),
                            alpha_prof_best=round(aq[0][0], 1) if aq else "", alpha_prof_n40=sum(1 for s_, _ in aq if s_ >= 40),
                            alpha_prof_where=where(aq), hmg_pfam_n=hp_n,
                            hmg_prof_best=round(hq_[0][0], 1) if hq_ else "", hmg_prof_n40=sum(1 for s_, _ in hq_ if s_ >= 40),
                            hmg_prof_where=where(hq_), secs=round(time.time() - t1, 1)))
            fh.flush()
            if i % 10 == 0:
                print(i, g, f"{time.time() - t0:.0f}s", file=sys.stderr, flush=True)


main()
