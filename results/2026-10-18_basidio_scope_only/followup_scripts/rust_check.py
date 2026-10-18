"""Independent bE/bW check on one rust genome: tblastn, miniprot and Pfam homeodomain, against the 4 rustHD records' proteins.

usage: rust_check.py ACC TAXID FASTA.gz OUT.json
"""
import json, os, re, subprocess, sys, tempfile, gzip, shutil
from collections import defaultdict
from pathlib import Path

acc, taxid, fa, out = sys.argv[1:5]
W = Path("/bigdata/stajichlab/jstajich/projects/mat-fam18")
refs = []  # (id, gene, seq)
for faa in sorted((W / "db/Basidiomycota/Pucciniales").glob("*/proteins.faa")):
    name = None; seq = []
    for line in faa.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if name: refs.append((name[0], name[1], "".join(seq)))
            m = re.match(r">(\S+?)\|gene_index=\d+\|name=(\w+)", line)
            name = (m.group(1).split("_rustHD")[0], m.group(2)) if m else None; seq = []
        else:
            seq.append(line.strip())
qlen = {f"{r[0]}|{r[1]}": len(r[2]) for r in refs}
tmp = Path(tempfile.mkdtemp(dir=os.environ.get("SCRATCH", "/tmp"), prefix="rc_"))
try:
    g = tmp / "g.fa"
    with gzip.open(fa, "rb") as fi, open(g, "wb") as fo: shutil.copyfileobj(fi, fo)
    q = tmp / "q.faa"
    q.write_text("".join(f">{r[0]}|{r[1]}\n{r[2]}\n" for r in refs))
    # contig lengths
    clen = {}; cur = None
    for line in open(g):
        if line.startswith(">"): cur = line[1:].split()[0]; clen[cur] = 0
        else: clen[cur] += len(line.strip())
    lens = sorted(clen.values(), reverse=True); tot = sum(lens); acc_ = 0; n50 = 0
    for l in lens:
        acc_ += l
        if acc_ >= tot / 2: n50 = l; break
    res = dict(acc=acc, taxid=taxid, genome_bp=tot, n_contigs=len(lens), n50=n50)
    subprocess.run(["makeblastdb", "-in", str(g), "-dbtype", "nucl", "-out", str(tmp / "db")], check=True, capture_output=True)
    tb = subprocess.run(["tblastn", "-query", str(q), "-db", str(tmp / "db"), "-evalue", "1e-3", "-num_threads", "2", "-max_target_seqs", "500",
                         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"], capture_output=True, text=True, check=True).stdout
    hsps = []
    for l in tb.splitlines():
        a = l.split("\t"); hsps.append(dict(q=a[0], gene=a[0].split("|")[1], contig=a[1], pid=float(a[2]), alen=int(a[3]), qs=int(a[4]), qe=int(a[5]),
                                            ss=int(a[6]), se=int(a[7]), ev=float(a[8]), bits=float(a[9])))
    # cluster HSPs per gene into loci (same contig, within 3 kb), keep the best-scoring loci
    loci = defaultdict(list)
    for gene in ("bE", "bW"):
        hs = sorted((h for h in hsps if h["gene"] == gene), key=lambda h: (h["contig"], min(h["ss"], h["se"])))
        cl = []
        for h in hs:
            lo, hi = min(h["ss"], h["se"]), max(h["ss"], h["se"])
            if cl and cl[-1]["contig"] == h["contig"] and lo - cl[-1]["hi"] <= 3000:
                cl[-1]["hi"] = max(cl[-1]["hi"], hi); cl[-1]["hs"].append(h)
            else:
                cl.append(dict(contig=h["contig"], lo=lo, hi=hi, hs=[h]))
        out_l = []
        for c in cl:
            best_q = defaultdict(lambda: [0, 0, 0.0]); cov = {}
            for h in c["hs"]:
                cov.setdefault(h["q"], []).append((h["qs"], h["qe"]))
            qc = 0; bq = None
            for qq, iv in cov.items():
                iv.sort(); m = 0; e = 0
                for s_, e_ in iv:
                    s_ = max(s_, e + 1); m += max(0, e_ - s_ + 1); e = max(e, e_)
                if m / qlen[qq] > qc: qc, bq = m / qlen[qq], qq
            bh = max(c["hs"], key=lambda h: h["bits"])
            strand = "+" if bh["se"] > bh["ss"] else "-"
            out_l.append(dict(contig=c["contig"], lo=c["lo"], hi=c["hi"], n_hsp=len(c["hs"]), best_bits=bh["bits"], best_pid=bh["pid"], qcov=round(qc, 2), strand=strand, best_query=bq))
        out_l.sort(key=lambda x: -x["best_bits"]); loci[gene] = out_l[:5]
    res["tblastn"] = dict(loci)
    # pair: a bE locus and a bW locus on one contig within 20 kb
    pairs = []
    for e in loci["bE"]:
        for w in loci["bW"]:
            if e["contig"] == w["contig"]:
                gap = max(w["lo"] - e["hi"], e["lo"] - w["hi"])
                if gap <= 20000:
                    left, right = (e, w) if e["lo"] <= w["lo"] else (w, e)
                    orient = "divergent" if (left["strand"] == "-" and right["strand"] == "+") else "convergent" if (left["strand"] == "+" and right["strand"] == "-") else "tandem"
                    pairs.append(dict(contig=e["contig"], gap=max(gap, 0), orient=orient, bE_bits=e["best_bits"], bW_bits=w["best_bits"], bE_qcov=e["qcov"], bW_qcov=w["qcov"],
                                      contig_len=clen[e["contig"]], bE_lo=e["lo"], bW_lo=w["lo"]))
    pairs.sort(key=lambda p: -(p["bE_bits"] + p["bW_bits"])); res["tblastn_pairs"] = pairs[:5]
    # miniprot
    mp = subprocess.run(["miniprot", "-t2", "--gff", str(g), str(q)], capture_output=True, text=True).stdout
    models = []
    cur = None
    for l in mp.splitlines():
        if l.startswith("#"): continue
        a = l.split("\t")
        if len(a) < 9: continue
        if a[2] == "mRNA":
            at = dict(x.split("=", 1) for x in a[8].split(";") if "=" in x)
            tq = at.get("Target", "").split()
            qq = tq[0]; qs, qe = int(tq[1]), int(tq[2])
            cur = dict(contig=a[0], start=int(a[3]), end=int(a[4]), strand=a[6], gene=qq.split("|")[1], q=qq, rank=int(at.get("Rank", 0)),
                       ident=float(at.get("Identity", 0)), pos=float(at.get("Positive", 0)), qcov=round((qe - qs + 1) / qlen[qq], 2), frameshift=0, stop=0, cds=0)
            models.append(cur)
        elif cur and a[2] == "CDS":
            cur["cds"] += 1
            if "Frameshift" in a[8]: cur["frameshift"] += 1
            if "Stop" in a[8] or "STOP" in a[8]: cur["stop"] += 1
        elif cur and a[2] == "stop_codon":
            pass
    for m in models:
        m["bp"] = m["end"] - m["start"] + 1
    best = {}
    for m in sorted(models, key=lambda m: -(m["ident"] * m["qcov"])):
        best.setdefault(m["gene"], []).append(m)
    res["miniprot_best"] = {k: v[:3] for k, v in best.items()}
    mp_pairs = []
    for e in best.get("bE", [])[:3]:
        for w in best.get("bW", [])[:3]:
            if e["contig"] == w["contig"]:
                gap = max(w["start"] - e["end"], e["start"] - w["end"])
                left, right = (e, w) if e["start"] <= w["start"] else (w, e)
                orient = "divergent" if (left["strand"] == "-" and right["strand"] == "+") else "convergent" if (left["strand"] == "+" and right["strand"] == "-") else "tandem"
                if gap <= 20000: mp_pairs.append(dict(contig=e["contig"], gap=max(gap, 0), orient=orient, bE_ident=e["ident"], bE_qcov=e["qcov"], bW_ident=w["ident"], bW_qcov=w["qcov"], contig_len=clen[e["contig"]]))
    res["miniprot_pairs"] = mp_pairs[:5]
    # Pfam homeodomain on six-frame ORFs of windows around the top tblastn loci
    import pyhmmer
    from Bio.Seq import Seq
    seqs = {}
    cur = None; buf = []
    want = {l["contig"] for gene in loci for l in loci[gene]}
    for line in open(g):
        if line.startswith(">"):
            if cur in want: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if cur in want: seqs[cur] = "".join(buf)
    abc = pyhmmer.easel.Alphabet.amino()
    nm = lambda x: x.decode() if isinstance(x, bytes) else x
    def mkseq(name, o, abc):
        try:
            return pyhmmer.easel.TextSequence(name=name.encode(), sequence=o).digitize(abc)
        except TypeError:
            return pyhmmer.easel.TextSequence(name=name, sequence=o).digitize(abc)
    hmms = []
    with pyhmmer.plan7.HMMFile("/bigdata/stajichlab/shared/lib/funannotate_db/Pfam-A.hmm") as hf:
        for h in hf:
            if nm(h.name) in ("Homeodomain", "Homeobox", "Homeobox_KN", "Homeodomain_2"): hmms.append(h)
            if len(hmms) == 4: break
    hd = []
    for gene in loci:
        for l in loci[gene][:3]:
            s = seqs[l["contig"]]; lo = max(0, l["lo"] - 4000); hi = min(len(s), l["hi"] + 4000)
            win = s[lo:hi]
            for strand, w_ in (("+", win), ("-", str(Seq(win).reverse_complement()))):
                for f in range(3):
                    prot = str(Seq(w_[f: len(w_) - ((len(w_) - f) % 3)]).translate())
                    for i, orf in enumerate(prot.split("*")):
                        if len(orf) < 40: continue
                        pass
                    # translate whole frame with stops as X so hmmsearch sees ORF fragments
                    frag = [o for o in prot.split("*") if len(o) >= 40]
                    seqblocks = [mkseq(f"{l['contig']}:{strand}{f}:{k}", o, abc) for k, o in enumerate(frag)]
                    if not seqblocks: continue
                    for hits in pyhmmer.hmmsearch(hmms, seqblocks, E=1e-3):
                        for hit in hits:
                            hd.append(dict(gene_locus=gene, contig=l["contig"], domain=nm(hits.query.name), score=round(hit.score, 1), evalue=hit.evalue))
    seen = {}
    for x in hd:
        k = (x["gene_locus"], x["contig"], x["domain"])
        if k not in seen or x["score"] > seen[k]["score"]: seen[k] = x
    res["pfam_homeodomain"] = list(seen.values())
    Path(out).write_text(json.dumps(res, indent=1))
finally:
    shutil.rmtree(tmp, ignore_errors=True)
