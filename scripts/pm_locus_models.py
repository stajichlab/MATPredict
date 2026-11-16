#!/usr/bin/env python3
"""Gene models and cross-genome orthology for HMG + homeodomain pairs (the PM-like locus).

Read-only discovery tool.

models : for each locus (TSV: name, genome fasta, contig, hd_start, hd_end), run exonerate protein2genome with the curated
         matPi proteins in a window around the homeodomain-like region, keep the best-scoring model, build the CDS from the
         exonerate `cds` features, translate it and run Pfam (HMG_box, HMG_box_2, Homeodomain, GA cutoff).
         Output: models.tsv (name, contig, strand, cds_start, cds_end, exons, aa_len, score, query, domains) and models.faa.
orthology : from a FASTA of protein sets (id = genome|role|gene, role in HMG, HD, HMG_locus, HD_locus), diamond all-vs-all;
         for each locus protein and each other genome, the best hit among that genome's proteins of the same domain class and
         whether it is that genome's locus protein (reciprocal best hit).
Usage: pm_locus_models.py models --loci LOCI.tsv --queries Q.faa --hmm HMM --out OUT
       pm_locus_models.py orthology --proteins ALL.faa --out OUT
"""
import argparse
import collections
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from region_pfam_check import read_fasta_region, revcomp, translate  # noqa: E402


def parse_exonerate_gff(text):
    """-> list of models; model = dict(score, query, strand, cds=[(start,end)], contig). Coordinates 1-based on the window."""
    models, cur = [], None
    for line in text.splitlines():
        if line.startswith("# --- START OF GFF DUMP"):
            cur = {"cds": [], "score": 0, "query": "", "strand": "+", "contig": ""}
        elif line.startswith("# --- END OF GFF DUMP"):
            if cur and cur["cds"]:
                models.append(cur)
            cur = None
        elif cur is not None and not line.startswith("#"):
            f = line.split("\t")
            if len(f) < 9:
                continue
            if f[2] == "gene":
                cur["score"] = float(f[5]) if f[5] != "." else 0
                cur["strand"] = f[6]
                cur["contig"] = f[0]
                for kv in f[8].split(";"):
                    kv = kv.strip()
                    if kv.startswith("sequence "):
                        cur["query"] = kv.split(" ", 1)[1]
            elif f[2] == "cds":
                cur["cds"].append((int(f[3]), int(f[4])))
    return models


def cds_protein(window_seq, cds, strand):
    """Translate the concatenated CDS exons of a model. cds: list of (start,end) 1-based on window_seq."""
    parts = [window_seq[s - 1:e] for s, e in sorted(cds)]
    dna = "".join(parts)
    if strand == "-":
        dna = revcomp(dna)
    prot = translate(dna)
    return prot.rstrip("*")


def do_models(a):
    os.makedirs(a.out, exist_ok=True)
    fo = open(os.path.join(a.out, "models.tsv"), "w")
    fa = open(os.path.join(a.out, "models.faa"), "w")
    fo.write("name\tcontig\tstrand\tcds_start\tcds_end\texons\taa_len\tscore\tquery\tdomains\n")
    for line in open(a.loci):
        name, fasta, contig, s, e = line.rstrip("\n").split("\t")
        ws, we = max(1, int(s) - a.pad), int(e) + a.pad
        seq = read_fasta_region(fasta, contig, ws, we)
        win = os.path.join(a.tmp, "win.fa")
        open(win, "w").write(f">{contig}\n{seq}\n")
        out = subprocess.run(["exonerate", "--model", "protein2genome", "--query", a.queries, "--target", win,
                              "--showtargetgff", "yes", "--showalignment", "no", "--showvulgar", "no", "--bestn", "1",
                              "--percent", str(a.percent), "--exhaustive", "no"], capture_output=True, text=True).stdout
        models = sorted(parse_exonerate_gff(out), key=lambda m: -m["score"])
        if not models:
            fo.write(f"{name}\t{contig}\t\t\t\t\t\t\t\tno_model\n")
            continue
        m = models[0]
        prot = cds_protein(seq, m["cds"], m["strand"])
        faa = os.path.join(a.tmp, "m.faa")
        open(faa, "w").write(f">m\n{prot}\n")
        tbl = os.path.join(a.tmp, "m.tbl")
        subprocess.run(["hmmsearch", "--cut_ga", "--domtblout", tbl, "-o", "/dev/null", a.hmm, faa], check=True)
        dom = sorted({l.split()[3] for l in open(tbl) if l[0] != "#"})
        starts = ws - 1 + min(x for x, _ in m["cds"])
        ends = ws - 1 + max(y for _, y in m["cds"])
        fo.write(f"{name}\t{contig}\t{m['strand']}\t{starts}\t{ends}\t{len(m['cds'])}\t{len(prot)}\t{m['score']}\t{m['query']}\t{','.join(dom) or 'none'}\n")
        fa.write(f">{name}\n{prot}\n")
    fo.close()
    fa.close()
    print(open(os.path.join(a.out, "models.tsv")).read())


def do_orthology(a):
    os.makedirs(a.out, exist_ok=True)
    seqs, cur = {}, None
    for line in open(a.proteins):
        if line[0] == ">":
            cur = line[1:].strip().split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(line.strip())
    meta = {k: k.split("|") for k in seqs}  # genome, role, gene
    db = os.path.join(a.tmp, "orth")
    subprocess.run(["diamond", "makedb", "--in", a.proteins, "-d", db, "--quiet"], check=True)
    tsv = os.path.join(a.tmp, "orth.tsv")
    subprocess.run(["diamond", "blastp", "-q", a.proteins, "-d", db, "-o", tsv, "--ultra-sensitive", "-e", "1e-3", "--max-target-seqs", "500",
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "--quiet"], check=True)
    hits = collections.defaultdict(list)
    for l in open(tsv):
        q, s, p, n, ev, b = l.rstrip("\n").split("\t")
        hits[q].append((float(b), s, float(p), int(n)))
    cls = {"HMG_locus": "HMG", "HD_locus": "HD", "HMG": "HMG", "HD": "HD"}
    fo = open(os.path.join(a.out, "orthology.tsv"), "w")
    fo.write("locus_protein\tother_genome\tbest_hit\thit_is_locus_protein\tidentity\taligned_aa\treciprocal\n")
    best = {}
    for q in seqs:
        for other in {m[0] for m in meta.values()}:
            cand = [(b, s, p, n) for b, s, p, n in hits[q] if meta[s][0] == other and cls[meta[s][1]] == cls[meta[q][1]] and s != q]
            best[(q, other)] = max(cand) if cand else None
    for q, (g, role, gene) in meta.items():
        if not role.endswith("_locus"):
            continue
        for other in sorted({m[0] for m in meta.values()} - {g}):
            b = best[(q, other)]
            if not b:
                fo.write(f"{q}\t{other}\tnone\t\t\t\t\n")
                continue
            s = b[1]
            back = best[(s, g)]
            fo.write(f"{q}\t{other}\t{s}\t{meta[s][1].endswith('_locus')}\t{b[2]}\t{b[3]}\t{bool(back and back[1] == q)}\n")
    fo.close()
    print(open(os.path.join(a.out, "orthology.tsv")).read())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["models", "orthology"])
    ap.add_argument("--loci")
    ap.add_argument("--queries")
    ap.add_argument("--hmm")
    ap.add_argument("--proteins")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pad", type=int, default=2500)
    ap.add_argument("--percent", type=int, default=10)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()
    do_models(a) if a.cmd == "models" else do_orthology(a)


if __name__ == "__main__":
    sys.exit(main())
