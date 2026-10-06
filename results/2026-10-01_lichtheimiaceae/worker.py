"""Per-genome Lichtheimiaceae MAT exploration worker (read-only analysis).

For each genome: miniprot the Mucoromycota reference proteins (sexM, sexP,
tptA, rnhA, glrA, algA, btbA) against the assembly, translate every alignment,
and score sexM/sexP-query models and (for LCG) all annotated proteins with the
shipped classifier HMMs (sexM, sexP, P1). Writes <out>/<genome_id>/{mp.tsv,
clf_models.tsv, clf_annot.tsv}.

Usage: worker.py genomes.tsv out_dir task_index n_tasks threads
"""
import csv, gzip, os, re, shutil, subprocess, sys
from pathlib import Path
import pyhmmer

E = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
HERE = Path(__file__).resolve().parent
manifest, out_dir, idx, n, threads = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
rows = list(csv.DictReader(open(manifest), delimiter="\t"))
mine = [r for i, r in enumerate(rows) if i % n == idx]
scratch = Path(os.environ.get("SCRATCH") or os.environ.get("TMPDIR") or "/tmp") / f"licht_{idx}"
scratch.mkdir(parents=True, exist_ok=True)

alpha = pyhmmer.easel.Alphabet.amino()
hmms = {}
for name in ("sexM", "sexP", "P1"):
    with pyhmmer.plan7.HMMFile(str(HERE / f"{name}.hmm")) as fh:
        hmms[name] = fh.read()


def score(seqs):
    """seqs: list of (id, protein). Returns {id: {sexM, sexP, P1}} best domain-independent seq scores."""
    if not seqs:
        return {}
    blk = pyhmmer.easel.DigitalSequenceBlock(alpha, [
        pyhmmer.easel.TextSequence(name=i.encode(), sequence=s.replace("*", "")).digitize(alpha)
        for i, s in seqs if len(s.replace("*", "")) >= 20])
    res = {i: {"sexM": 0.0, "sexP": 0.0, "P1": 0.0} for i, _ in seqs}
    for name, hmm in hmms.items():
        for top in pyhmmer.hmmsearch([hmm], blk, cpus=int(threads), T=0.0, domT=0.0):
            for hit in top:
                nm = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                res[nm][name] = max(res[nm][name], hit.score)
    return res



CODON = {}
_b = "TCAG"; _aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_b):
    for _j, _c in enumerate(_b):
        for _k, _d in enumerate(_b):
            CODON[_a + _c + _d] = _aa[16 * _i + 4 * _j + _k]


def translate(nt):
    nt = nt.upper()
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def read_fasta(path):
    seqs, nm, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if nm: seqs[nm] = "".join(buf)
            nm = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if nm: seqs[nm] = "".join(buf)
    return seqs


def tblastn_exonerate(genome, gid, od):
    db = scratch / f"{gid}.db"
    subprocess.run([f"{E}/makeblastdb", "-in", genome, "-dbtype", "nucl", "-out", str(db)], capture_output=True, check=True)
    t = subprocess.run([f"{E}/tblastn", "-query", str(HERE / "sexMP_refs.faa"), "-db", str(db), "-evalue", "1e-3",
                        "-seg", "no", "-num_threads", threads,
                        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"],
                       capture_output=True, text=True)
    hits = []
    for line in t.stdout.splitlines():
        q, s, pid, ln, qs, qe, ss, se, ev, bs = line.split("\t")
        a, b = sorted((int(ss), int(se)))
        hits.append((s, a, b, q, float(bs), float(pid)))
    hits.sort()
    loci = []
    for h in hits:
        if loci and loci[-1]["contig"] == h[0] and h[1] <= loci[-1]["end"] + 3000:
            L = loci[-1]; L["end"] = max(L["end"], h[2]); L["hits"].append(h)
        else:
            loci.append(dict(contig=h[0], start=h[1], end=h[2], hits=[h]))
    gseq = None
    refs = read_fasta(HERE / "sexMP_refs.faa")
    out = []
    rows = []
    for li, L in enumerate(loci):
        best = {}
        for h in L["hits"]:
            best[h[3]] = max(best.get(h[3], 0), h[4])
        top = sorted(best, key=lambda q: -best[q])[:3]  # top refs for this locus
        if gseq is None:
            gseq = read_fasta(genome)
        cs = gseq[L["contig"]]
        a = max(0, L["start"] - 5000); b = min(len(cs), L["end"] + 5000)
        reg = scratch / f"{gid}_reg.fa"; reg.write_text(f">r\n{cs[a:b]}\n")
        bestm = None
        for q in top:
            qf = scratch / f"{gid}_q.fa"; qf.write_text(f">q\n{refs[q]}\n")
            ex = subprocess.run([f"{E}/exonerate", "-m", "protein2genome", "--showalignment", "no", "--showvulgar", "no",
                                 "--bestn", "1", "--ryo", "##R\t%s\t%tab\t%tae\t%tS\n%tcs\n", str(qf), str(reg)],
                                capture_output=True, text=True, timeout=300)
            txt = ex.stdout
            m = re.search(r"##R\t(\d+)\t(\d+)\t(\d+)\t([+-])\n([ACGTNacgtn\n]+)", txt)
            if not m: continue
            sc = int(m.group(1)); cds = m.group(5).replace("\n", "")
            if bestm is None or sc > bestm[0]:
                bestm = (sc, q, a + int(m.group(2)), a + int(m.group(3)), m.group(4), translate(cds))
        mid = f"{gid}|L{li}"
        tb = max(best.values())
        if bestm:
            out.append((mid, bestm[5]))
        rows.append(dict(mid=mid, contig=L["contig"], start=L["start"], end=L["end"], top_ref=top[0], top_tblastn_bits=tb,
                         n_hits=len(L["hits"]), model_ref=bestm[1] if bestm else "", model_start=bestm[2] if bestm else "",
                         model_end=bestm[3] if bestm else "", model_strand=bestm[4] if bestm else "",
                         model_len=len(bestm[5]) if bestm else 0, exonerate_score=bestm[0] if bestm else ""))
    with open(od / "loci.tsv", "w") as fo:
        cols = ["mid", "contig", "start", "end", "top_ref", "top_tblastn_bits", "n_hits", "model_ref", "model_start", "model_end", "model_strand", "model_len", "exonerate_score"]
        fo.write("\t".join(cols) + "\n")
        for r_ in rows:
            fo.write("\t".join(str(r_[c]) for c in cols) + "\n")
    with open(od / "models.faa", "w") as fo:
        for i, sq in out:
            fo.write(f">{i}\n{sq}\n")
    for f in scratch.glob(f"{gid}.db*"):
        f.unlink()
    return out


for r in mine:
    gid = r["genome_id"]
    od = out_dir / gid
    if (od / "done").exists():
        continue
    od.mkdir(parents=True, exist_ok=True)
    g = r["genome"]
    if g.endswith(".gz"):
        local = scratch / f"{gid}.fa"
        with gzip.open(g, "rb") as fi, open(local, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        g = str(local)
    # miniprot: gff + translated sequences
    p = subprocess.run([f"{E}/miniprot", "-t", threads, "--gff", "--trans", "-p", "0.5", "-N", "20", g, str(HERE / "refs.faa")],
                       capture_output=True, text=True)
    models = []  # (mid, query, gene, contig, start, end, strand, identity, positive, protein)
    cur = None
    mid = 0
    for line in p.stdout.splitlines():
        if line.startswith("##PAF"):
            f = line.split("\t")
            q = f[1]; gene = q.split("|")[-1]
            tags = {t.split(":", 2)[0]: t.split(":", 2)[2] for t in f[13:] if t.count(":") >= 2}
            nmatch, alen = int(f[10]), int(f[11]) if f[11].isdigit() else 0
            cur = dict(query=q, gene=gene, contig=f[6], start=int(f[8]), end=int(f[9]), strand=f[5],
                       qlen=int(f[2]), qcov=(int(f[4]) - int(f[3])) / max(1, int(f[2])),
                       identity=float(nmatch) / max(1, alen))
        elif line.startswith("##STA") and cur is not None:
            mid += 1
            cur["mid"] = f"{gid}|m{mid}"
            cur["protein"] = line.split("\t", 1)[1].strip()
            models.append(cur)
            cur = None
    with open(od / "mp.tsv", "w") as fo:
        cols = ["mid", "query", "gene", "contig", "start", "end", "strand", "qlen", "qcov", "identity", "plen"]
        fo.write("\t".join(cols) + "\n")
        for m in models:
            fo.write("\t".join(str(m.get(c, len(m["protein"]) if c == "plen" else "")) for c in cols) + "\n")
    # sexM/sexP: miniprot misses these divergent HMG genes, so use tblastn loci + exonerate models
    core = tblastn_exonerate(g, gid, od)
    sc = score(core)
    with open(od / "clf_models.tsv", "w") as fo:
        fo.write("mid\tsexM\tsexP\tP1\n")
        for i, _ in core:
            s = sc.get(i, {})
            fo.write(f"{i}\t{s.get('sexM',0):.1f}\t{s.get('sexP',0):.1f}\t{s.get('P1',0):.1f}\n")
    if r["proteins"]:
        seqs = []
        nm = None; buf = []
        for line in open(r["proteins"]):
            if line.startswith(">"):
                if nm: seqs.append((nm, "".join(buf)))
                nm = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if nm: seqs.append((nm, "".join(buf)))
        sc = score(seqs)
        with open(od / "clf_annot.tsv", "w") as fo:
            fo.write("protein\tsexM\tsexP\tP1\n")
            for i, s in sc.items():
                if max(s.values()) >= 20:
                    fo.write(f"{i}\t{s['sexM']:.1f}\t{s['sexP']:.1f}\t{s['P1']:.1f}\n")
    if g.startswith(str(scratch)):
        os.remove(g)
    (od / "done").write_text("ok\n")
    print("done", gid, len(models), flush=True)
