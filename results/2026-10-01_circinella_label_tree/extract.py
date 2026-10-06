"""Extract the sex-locus HMG protein for the Circinella group and reference genomes.

For each genome: take the best miniprot rnhA hit (from the Lichtheimiaceae
exploration's mp.tsv), find sexM/sexP miniprot hits on the same contig within
15 kb of rnhA (the Circinella-group arrangement). If none, take the best
sexM/sexP miniprot hit anywhere (reference genomes with other layouts).
Re-run miniprot of all 18 curated sexM/sexP proteins on the locus region
(+-8 kb) with --trans and keep the best-scoring full model. For LCG genomes,
also take the annotated protein overlapping the locus (funannotate gff3).
Score both with the shipped classifier (sexM, sexP, P1).

Read-only. Usage: extract.py genomes.tsv out_prefix
"""
import csv, gzip, re, subprocess, sys
from pathlib import Path
import pyhmmer

E = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
R = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/results")
LICHT = R / "2026-10-01_lichtheimiaceae"
CLF = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db/Mucoromycota/classifiers/MAT")
REFS = LICHT / "sexMP_refs.faa"
SCR = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("/tmp")

alpha = pyhmmer.easel.Alphabet.amino()
hmms = {}
for name, p in (("sexM", CLF / "sexM.hmm"), ("sexP", CLF / "sexP.hmm"), ("P1", CLF / "paralogs" / "P1.hmm")):
    with pyhmmer.plan7.HMMFile(str(p)) as fh:
        hmms[name] = fh.read()


def score(seq):
    s = seq.replace("*", "").replace("X", "")
    if len(s) < 20:
        return {"sexM": 0.0, "sexP": 0.0, "P1": 0.0}
    blk = pyhmmer.easel.DigitalSequenceBlock(alpha, [pyhmmer.easel.TextSequence(name=b"q", sequence=s).digitize(alpha)])
    out = {}
    for n, h in hmms.items():
        best = 0.0
        for top in pyhmmer.hmmsearch([h], blk, T=0.0, domT=0.0):
            for hit in top:
                best = max(best, hit.score)
        out[n] = round(best, 1)
    return out


def read_contig(path, contig):
    op = gzip.open if str(path).endswith(".gz") else open
    seq, on = [], False
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on:
                    break
                on = line[1:].split()[0] == contig
            elif on:
                seq.append(line.strip())
    return "".join(seq)


def read_fasta(path):
    seqs, nm, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if nm:
                seqs[nm] = "".join(buf)
            nm = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if nm:
        seqs[nm] = "".join(buf)
    return seqs


def mp_rows(gid):
    p = LICHT / "out" / gid / "mp.tsv"
    if not p.exists():
        return []
    return list(csv.DictReader(open(p), delimiter="\t"))


def annotated(gff, prots, contig, a, b):
    """Protein of the mRNA overlapping [a,b] most."""
    if not gff or not Path(gff).exists():
        return None, None
    best = (0, None)
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != contig or f[2] != "mRNA":
            continue
        s, e = int(f[3]), int(f[4])
        ov = min(e, b) - max(s, a)
        if ov > best[0]:
            m = re.search(r"ID=([^;]+)", f[8])
            best = (ov, m.group(1) if m else None)
    if not best[1]:
        return None, None
    ps = read_fasta(prots) if prots and Path(prots).exists() else {}
    return best[1], ps.get(best[1])


def miniprot_locus(genome, gid, contig, a, b):
    cs = read_contig(genome, contig)
    if not cs:
        return None
    lo, hi = max(0, a - 8000), min(len(cs), b + 8000)
    reg = SCR / f"{gid}_reg.fa"
    reg.write_text(f">{contig}_{lo}\n{cs[lo:hi]}\n")
    p = subprocess.run([f"{E}/miniprot", "--trans", "-p", "0.3", "-N", "5", str(reg), str(REFS)],
                       capture_output=True, text=True)
    best, cur = None, None
    for line in p.stdout.splitlines():
        if not line.startswith("#") and "\t" in line:
            f = ["PAF"] + line.split("\t")
            ms = [x for x in f if x.startswith("AS:i:")]
            cur = {"query": f[1], "qlen": int(f[2]), "qs": int(f[3]), "qe": int(f[4]),
                   "ts": lo + int(f[8]), "te": lo + int(f[9]), "AS": int(ms[0][5:]) if ms else 0}
        elif line.startswith("##STA") and cur:
            cur["protein"] = line.split("\t")[1].strip()
            if best is None or cur["AS"] > best["AS"]:
                best = cur
            cur = None
    return best


rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
out = Path(sys.argv[2])
fa = open(str(out) + ".faa", "w")
tab = open(str(out) + ".tsv", "w")
cols = ["gid", "group", "species", "label", "rnhA_contig", "rnhA_start", "rnhA_end", "locus_contig", "locus_start",
        "locus_end", "dist_to_rnhA_kb", "near_rnhA", "mp_query", "mp_len", "mp_qcov", "mp_sexM", "mp_sexP", "mp_P1",
        "annot_id", "annot_len", "annot_sexM", "annot_sexP", "annot_P1", "chosen", "chosen_len"]
tab.write("\t".join(cols) + "\n")
for r in rows:
    gid = r["gid"]
    mps = mp_rows(gid)
    rn = [m for m in mps if m["gene"] == "rnhA"]
    rbest = max(rn, key=lambda m: float(m["qcov"]) * float(m["identity"])) if rn else None
    sx = [m for m in mps if m["gene"] in ("sexM", "sexP")]
    near = []
    if rbest:
        rs, re_ = int(rbest["start"]), int(rbest["end"])
        for m in sx:
            if m["contig"] == rbest["contig"]:
                d = max(0, max(int(m["start"]), rs) - min(int(m["end"]), re_))
                if d <= 15000:
                    near.append((d, m))
    if near:
        d, loc = min(near, key=lambda t: (t[0], -float(t[1]["qcov"]) * float(t[1]["identity"])))
        nr = True
    elif sx:
        loc = max(sx, key=lambda m: float(m["qcov"]) * float(m["identity"])); d = None; nr = False
    else:
        tab.write("\t".join([gid, r["group"], r["species"], r["label"]] + [""] * (len(cols) - 4)) + "\n")
        continue
    c, a, b = loc["contig"], int(loc["start"]), int(loc["end"])
    mp = miniprot_locus(r["genome"], gid, c, a, b)
    aid, aseq = annotated(r.get("gff"), r.get("proteins"), c, a, b)
    ms = score(mp["protein"]) if mp and mp.get("protein") else {}
    asc = score(aseq) if aseq else {}
    # choose: annotated if it carries an HMG signal and is at least as long; else miniprot
    choice, cseq = None, None
    mlen = len(mp.get("protein", "")) if mp else 0
    if aseq and max(asc.get("sexM", 0), asc.get("sexP", 0)) >= 20 and (not mp or 0.9 * mlen <= len(aseq) <= 1.5 * mlen):
        choice, cseq = "annot", aseq.replace("*", "")
    elif mp and mp.get("protein"):
        choice, cseq = "miniprot", mp["protein"].replace("*", "")
    if cseq:
        fa.write(f">{gid}\n{cseq}\n")
    vals = [gid, r["group"], r["species"], r["label"],
            rbest["contig"] if rbest else "", rbest["start"] if rbest else "", rbest["end"] if rbest else "",
            c, a, b, "" if d is None else round(d / 1000, 1), nr,
            mp["query"] if mp else "", len(mp["protein"]) if mp and mp.get("protein") else "",
            round((mp["qe"] - mp["qs"]) / mp["qlen"], 2) if mp else "",
            ms.get("sexM", ""), ms.get("sexP", ""), ms.get("P1", ""),
            aid or "", len(aseq) if aseq else "", asc.get("sexM", ""), asc.get("sexP", ""), asc.get("P1", ""),
            choice or "", len(cseq) if cseq else ""]
    tab.write("\t".join(str(v) for v in vals) + "\n")
    tab.flush(); fa.flush()
fa.close(); tab.close()
