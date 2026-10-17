#!/usr/bin/env python
"""Triage of Ascomycota genomes the detector did not call: is there a MAT-like gene, and is it linked to the flanks?

    triage_uncalled.py run --genome-id G --genomes-dir DIR --queries refs.faa --profiles MAT_core --out row.json
    triage_uncalled.py summarize ROW.json [...] [--groups groups.tsv]

For one genome: translate all six frames, split at stops (segments of at least MIN_AA residues) and search them with the
alpha-box (MAT1-1-1) and HMG-box (MAT1-2-1) profiles, which is more sensitive to divergent domains than a
protein-to-genome alignment; miniprot places the conserved SLA2, APN2 and COX13 flank genes from the curated proteins.
The profile search is independent of the pipeline's own search. Classes (independent of the pipeline's call):

  C1  no MAT-like gene: no MAT1-1-1 or MAT1-2-1 hit scores on its own profile (below WEAK bits)
  C2  MAT-like gene present, no flank gene (SLA2, APN2, COX13) within LINK_BP on the same contig
  C3  MAT-like gene linked to a flank gene: if the pipeline did not call it, a pipeline miss
  C0  no flank gene found either (fragmented or divergent flanks): cannot say

It also records the SLA2-APN2 interval when both flanks are on one contig (a conserved flank block with no gene
between them is the Xylariales-type "empty flank block"). Thresholds are provisional, taken from the curated records'
leave-one-genus-out scores, and every raw score is kept so they can be moved.
"""
import argparse
import csv
import gzip
import json
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

STRONG, WEAK = 100.0, 40.0      # profile bits: strong, and the floor under which a hit is not MAT-like
MIN_ID, MIN_COV = 0.30, 0.40    # miniprot identity and query coverage for a flank gene
MIN_AA = 30                     # shortest six-frame segment searched
LINK_BP = 30_000
CORE = {"MAT1-1-1": "MAT1-1", "MAT1-2-1": "MAT1-2"}
FLANKS = {"SLA2", "APN2", "COX13", "sla2", "apn2", "cox13"}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    seqs, cur = {}, None
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0]; seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def gene_of(header: str) -> str:
    return header.split("|")[-1]          # record|geneN|NAME


def parse_miniprot(gff: str):
    """mRNA records with their CDS exons: dict(target, ident, contig, strand, start, end, exons, qstart, qend)."""
    out, cur = [], None
    for line in gff.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] == "mRNA":
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            t = attrs.get("Target", "").split()
            cur = dict(target=t[0], qstart=int(t[1]), qend=int(t[2]), ident=float(attrs.get("Identity", 0)),
                       contig=f[0], strand=f[6], start=int(f[3]), end=int(f[4]), exons=[], id=attrs["ID"],
                       frameshift="Frameshift" in attrs, stop="StopCodon" in attrs)
            out.append(cur)
        elif f[2] == "CDS" and cur is not None:
            cur["exons"].append((int(f[3]), int(f[4])))
    return out


def translate(exons, strand, seq):
    from Bio.Seq import Seq
    dna = "".join(seq[a - 1:b] for a, b in sorted(exons))
    if strand == "-":
        dna = dna.translate(_COMP)[::-1]
    return str(Seq(dna[: len(dna) // 3 * 3]).translate()).replace("*", "X")


def load_profiles(d: Path):
    import pyhmmer
    out = {}
    for gene in CORE:
        with pyhmmer.plan7.HMMFile(str(d / f"{gene}.hmm")) as fh:
            out[gene] = fh.read()
    return out


def best_scores(profiles, protein: str) -> dict:
    from MATPredict.detect.classifier_build import score
    return {g: score(h, {"p": protein})["p"] for g, h in profiles.items()}


def open_genome(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        from Bio import SeqIO
        return {r.id: str(r.seq).upper() for r in SeqIO.parse(fh, "fasta")}


def sixframe_hits(genome: dict, profiles: dict, keep: int = 6) -> list[dict]:
    """Best segment per profile and contig-region: (gene, contig, start, end, strand, aa, alpha, hmg, own)."""
    import pyhmmer
    from Bio.Seq import Seq
    alphabet = pyhmmer.easel.Alphabet.amino()
    segs, meta = [], {}
    for contig, seq in genome.items():
        for strand in "+-":
            dna = seq if strand == "+" else seq.translate(_COMP)[::-1]
            for f in range(3):
                prot = str(Seq(dna[f:][: (len(dna) - f) // 3 * 3]).translate())
                pos = 0
                for piece in prot.split("*"):
                    if len(piece) >= MIN_AA:
                        a0, a1 = f + 3 * pos, f + 3 * (pos + len(piece))
                        name = f"s{len(segs)}"
                        if strand == "+":
                            lo, hi = a0 + 1, a1
                        else:
                            lo, hi = len(dna) - a1 + 1, len(dna) - a0
                        segs.append(pyhmmer.easel.TextSequence(name=name.encode(), sequence=piece.replace("X", "A")).digitize(alphabet))
                        meta[name] = (contig, lo, hi, strand, len(piece))
                    pos += len(piece) + 1
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, segs)
    pipe = pyhmmer.plan7.Pipeline(alphabet, Z=1, E=1e9, domE=1e9, bias_filter=False, F1=1.0, F2=1.0, F3=1.0)
    best: dict[str, dict] = {}
    for gene, hmm in profiles.items():
        for hit in pipe.search_hmm(hmm, block):
            name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            sc = best.setdefault(name, {"MAT1-1-1": 0.0, "MAT1-2-1": 0.0})
            sc[gene] = max(sc[gene], float(hit.score))
    rows = []
    for name, sc in best.items():
        contig, lo, hi, strand, aa = meta[name]
        gene = "MAT1-1-1" if sc["MAT1-1-1"] >= sc["MAT1-2-1"] else "MAT1-2-1"
        own = sc[gene]
        if own >= WEAK:
            rows.append(dict(gene=gene, contig=contig, start=lo, end=hi, strand=strand, aa=aa,
                             alpha=round(sc["MAT1-1-1"], 1), hmg=round(sc["MAT1-2-1"], 1), own=round(own, 1)))
    return sorted(rows, key=lambda r: -r["own"])[:keep * 4]


def run(a) -> int:
    gdir = Path(a.genomes_dir)
    path = next((p for p in (gdir / f"{a.genome_id}.fa.gz", gdir / f"{a.genome_id}.fa") if p.exists()), None)
    out = dict(genome=a.genome_id)
    if path is None:
        out["klass"] = "no_genome"; Path(a.out).write_text(json.dumps(out)); return 0
    genome = open_genome(path)
    out["n_contigs"], out["length_bp"] = len(genome), sum(map(len, genome.values()))
    queries = read_fasta(a.queries)
    flank_q = {k: v for k, v in queries.items() if gene_of(k) in FLANKS}
    qlen = {k: len(v) for k, v in flank_q.items()}
    profiles = load_profiles(Path(a.profiles))
    core_hits = sixframe_hits(genome, profiles)
    with tempfile.TemporaryDirectory() as tmp:
        fna, qfa = Path(tmp) / "g.fna", Path(tmp) / "flank.faa"
        with open(fna, "w") as fh:
            for k, v in genome.items():
                fh.write(f">{k}\n{v}\n")
        with open(qfa, "w") as fh:
            for k, v in flank_q.items():
                fh.write(f">{k}\n{v}\n")
        p = subprocess.run(["miniprot", "-t", str(a.threads), "--gff", "-p", "0.3", "-N", "5", str(fna), str(qfa)],
                           capture_output=True, text=True)
    if p.returncode != 0:
        out["klass"] = "miniprot_error"; out["returncode"] = p.returncode; out["stderr"] = p.stderr[-300:]
        Path(a.out).write_text(json.dumps(out)); return 0
    flank_hits = []
    for h in parse_miniprot(p.stdout):
        gene = gene_of(h["target"])
        cov = (h["qend"] - h["qstart"] + 1) / max(qlen.get(h["target"], 1), 1)
        if gene in FLANKS and h["ident"] >= MIN_ID and cov >= MIN_COV:
            flank_hits.append(dict(gene=gene.upper(), contig=h["contig"], start=h["start"], end=h["end"],
                                   ident=round(h["ident"], 3), cov=round(cov, 2)))
    for c in core_hits:
        c["mat_like"] = c["own"] >= WEAK
        c["linked"] = any(f["contig"] == c["contig"] and (f["start"] <= c["end"] and f["end"] >= c["start"]
                          or min(abs(f["start"] - c["end"]), abs(f["end"] - c["start"])) <= LINK_BP) for f in flank_hits)
    out["core_hits"] = sorted(core_hits, key=lambda c: -c["own"])[:8]
    out["flank_hits"] = flank_hits[:20]
    by_flank = defaultdict(list)
    for f in flank_hits:
        by_flank[f["gene"]].append(f)
    block = None
    for s in by_flank.get("SLA2", []):
        for ap in by_flank.get("APN2", []):
            if s["contig"] == ap["contig"]:
                lo, hi = sorted([s["end"], ap["start"]]) if s["end"] < ap["start"] else sorted([ap["end"], s["start"]])
                if block is None or abs(hi - lo) < block["gap_bp"]:
                    block = dict(contig=s["contig"], gap_bp=hi - lo)
    out["sla2_apn2_interval"] = block
    mat_like = [c for c in core_hits if c["mat_like"]]
    if mat_like and any(c["linked"] for c in mat_like):
        out["klass"] = "C3_linked"
    elif mat_like:
        out["klass"] = "C2_unlinked" if flank_hits else "C0_mat_like_no_flank"
    else:
        out["klass"] = "C1_no_mat_like" if flank_hits else "C0_no_flank"
    out["strong"] = any(c["own"] >= STRONG for c in core_hits)
    out["strong_linked"] = any(c["own"] >= STRONG and c["linked"] for c in core_hits)
    Path(a.out).write_text(json.dumps(out, indent=1))
    return 0


def summarize(a) -> int:
    rows = [json.loads(Path(p).read_text()) for p in a.rows]
    groups = {}
    if a.groups:
        groups = {r["genome"]: r for r in csv.DictReader(open(a.groups), delimiter="\t")}
    from collections import Counter
    tab = defaultdict(Counter)
    for r in rows:
        g = groups.get(r["genome"], {})
        tab[(g.get("group", "?"), g.get("order", "?"))][r["klass"]] += 1
        if r.get("strong_linked"):
            tab[(g.get("group", "?"), g.get("order", "?"))]["strong_linked"] += 1
    print(f"{'group':9s} {'order':20s} {'n':>4s}  C3_linked  (strong)  C2_unlinked  C1_no_mat_like  C0  other")
    for (grp, order), c in sorted(tab.items()):
        n = sum(v for k, v in c.items() if k != "strong_linked"); c0 = c["C0_no_flank"] + c["C0_mat_like_no_flank"]
        other = n - c["C3_linked"] - c["C2_unlinked"] - c["C1_no_mat_like"] - c0
        print(f"{grp:9s} {order:20s} {n:4d}  {c['C3_linked']:9d}  {c['strong_linked']:8d}  {c['C2_unlinked']:11d}  {c['C1_no_mat_like']:14d}  {c0:2d}  {other}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--genome-id", required=True); r.add_argument("--genomes-dir", required=True)
    r.add_argument("--queries", required=True); r.add_argument("--profiles", required=True)
    r.add_argument("--out", required=True); r.add_argument("--threads", type=int, default=2)
    s = sub.add_parser("summarize"); s.add_argument("rows", nargs="+"); s.add_argument("--groups")
    a = ap.parse_args(argv)
    return run(a) if a.cmd == "run" else summarize(a)


if __name__ == "__main__":
    sys.exit(main())
