"""Panel metric for helixer-lite window models against the pipeline and the NCBI annotation.

Pre-registered in docs/notes/2026-10-15_session-handoff.md (metric section). Two stages:

  helixer  run helixer-lite on the locus ±5 kb windows of one genome; write its genes as TSV
  score    compare, per genome, the pipeline genes, the helixer genes and the NCBI genes

Match: same contig and strand, overlap at least 50% of the longer of the two genes.
Boundary agreement: both ends within 30 bp. Strata: the pipeline gene's role (core MAT, flank, other);
NCBI genes the pipeline does not model are counted as 'unmodelled'.

usage:
  panel_helixer_metric.py helixer --genome G.fna --detected detected_loci.gff3 --out helixer.tsv \
      --helixerlite /path/helixerlite [--cpus 8]
  panel_helixer_metric.py score --run NAME detected_loci.gff3 helixer.tsv ncbi.gff [--run ...]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

TOL_BP = 30
MIN_RECIPROCAL = 0.5


def read_detected_genes(path: Path):
    """(contig, start, end, strand, name, role) for present genes, and the MAT_locus spans."""
    genes, spans = [], []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "MAT_locus":
            spans.append((f[0], int(f[3]), int(f[4])))
        elif f[2] == "gene" and attrs.get("present") == "true":
            genes.append((f[0], int(f[3]), int(f[4]), f[6], attrs.get("Name", "?"), attrs.get("role", "")))
    return genes, spans


def read_ncbi_genes(path: Path):
    out = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        out.append((f[0], int(f[3]), int(f[4]), f[6], attrs.get("Name", attrs.get("ID", "?"))))
    return out


def match(a, b) -> int:
    """Overlap in bp if a and b match under the 50%-reciprocal rule on the same strand, else 0."""
    if a[0] != b[0] or a[3] != b[3]:
        return 0
    ov = min(a[2], b[2]) - max(a[1], b[1]) + 1
    if ov <= 0:
        return 0
    longer = max(a[2] - a[1] + 1, b[2] - b[1] + 1)
    return ov if ov >= MIN_RECIPROCAL * longer else 0


def best(ref, cands):
    best_c, best_ov = None, 0
    for c in cands:
        ov = match(ref, c)
        if ov > best_ov:
            best_c, best_ov = c, ov
    return best_c


def within(a, b) -> bool:
    return abs(a[1] - b[1]) <= TOL_BP and abs(a[2] - b[2]) <= TOL_BP


def stratum(role: str) -> str:
    if role == "core_MAT":
        return "core"
    if role.startswith("flanking"):
        return "flank"
    return "other"


def run_helixer(args) -> None:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from Bio import SeqIO
    from MATPredict.detect.helixer_windows import run_helixer_windows

    _, spans = read_detected_genes(Path(args.detected))
    seqs, lengths = {}, {}
    wanted = {c for c, _, _ in spans}
    for rec in SeqIO.parse(args.genome, "fasta"):
        if rec.id in wanted:
            seqs[rec.id] = str(rec.seq).upper()
            lengths[rec.id] = len(rec.seq)
    track = run_helixer_windows(spans, seqs, lengths, flank=5000, lineage="fungi", cpus=args.cpus,
                                binary=[args.helixerlite])
    with open(args.out, "w") as fh:
        fh.write("contig\tstart\tend\tstrand\tname\n")
        for g in track.genes:
            fh.write(f"{g.contig}\t{g.start}\t{g.end}\t{g.strand}\t{g.name}\n")
    print(f"helixer: {len(track.genes)} models in {len(spans)} locus span(s) -> {args.out}")


def read_helixer(path: Path):
    out = []
    for line in path.read_text().splitlines()[1:]:
        f = line.split("\t")
        out.append((f[0], int(f[1]), int(f[2]), f[3], f[4]))
    return out


def score_one(name, detected, helixer_tsv, ncbi_gff):
    ours, spans = read_detected_genes(Path(detected))
    ours_t = [(g[0], g[1], g[2], g[3], g[4]) for g in ours]
    role_of = {(g[0], g[1], g[2], g[3]): g[5] for g in ours}
    hx = read_helixer(Path(helixer_tsv))
    ncbi = read_ncbi_genes(Path(ncbi_gff))

    counts = defaultdict(Counter)
    ours_in_ncbi = set()
    for n in ncbi:
        o = best(n, ours_t)
        if o is None:
            # NCBI gene the pipeline does not model: counted only when it lies in a locus window
            if any(c == n[0] and s <= n[2] and e >= n[1] for c, s, e in spans):
                counts["unmodelled"]["ncbi"] += 1
                if best(n, hx):
                    counts["unmodelled"]["helixer_matched"] += 1
            continue
        st = stratum(role_of[(o[0], o[1], o[2], o[3])])
        counts[st]["ncbi"] += 1
        ours_in_ncbi.add((o[0], o[1], o[2], o[3]))
        if within(n, o):
            counts[st]["pipeline_bounded"] += 1
        h = best(n, hx)
        if h is not None:
            counts[st]["helixer_matched"] += 1
            if within(n, h):
                counts[st]["helixer_bounded"] += 1
    unmatched_hx = sum(1 for h in hx if best(h, ncbi) is None)
    print(f"== {name}  (pipeline genes {len(ours_t)}, helixer models {len(hx)}, NCBI genes {len(ncbi)})")
    print(f"   {'stratum':10s} {'n_ncbi':>6s} {'pipe_bounded':>13s} {'hx_matched':>11s} {'hx_bounded':>11s}")
    for st in ("core", "flank", "other", "unmodelled"):
        c = counts[st]
        n = c["ncbi"]
        if n == 0 and not c["helixer_matched"]:
            continue
        pb = f"{c['pipeline_bounded']}/{n}" if st != "unmodelled" else "-"
        hm = f"{c['helixer_matched']}/{n}" if n else f"{c['helixer_matched']}"
        hb = f"{c['helixer_bounded']}/{c['helixer_matched']}" if c["helixer_matched"] and st != "unmodelled" else "-"
        print(f"   {st:10s} {n:6d} {pb:>13s} {hm:>11s} {hb:>11s}")
    print(f"   helixer models matching no NCBI gene (window extras): {unmatched_hx}")
    return counts


def run_score(args) -> None:
    for name, detected, helixer_tsv, ncbi_gff in args.run:
        score_one(name, detected, helixer_tsv, ncbi_gff)


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="stage", required=True)
    h = sub.add_parser("helixer")
    h.add_argument("--genome", required=True)
    h.add_argument("--detected", required=True)
    h.add_argument("--out", required=True)
    h.add_argument("--helixerlite", required=True)
    h.add_argument("--cpus", type=int, default=4)
    s = sub.add_parser("score")
    s.add_argument("--run", nargs=4, metavar=("NAME", "DETECTED", "HELIXER_TSV", "NCBI_GFF"),
                   action="append", required=True)
    args = p.parse_args(argv)
    if args.stage == "helixer":
        run_helixer(args)
    else:
        run_score(args)


if __name__ == "__main__":
    main()
