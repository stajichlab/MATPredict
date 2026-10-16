"""Score the flank rule with and without the query check, per genome, against the NCBI annotation.

For each flanking gene of a detection that matches an NCBI gene: baseline (homology ends), the rule without the
query check, and the rule with it (helixer's protein must align the curated query protein better than the
homology protein does). Prints one line per gene and the totals. Metric definitions are those of
panel_helixer_metric.py (50% reciprocal match, 30 bp tolerance).

usage:
  panel_rule_check.py --db-root DB --run NAME GENOME.fna DETECT_DIR NCBI_GFF HELIXER_GFF [--run ...] --helixerlite EXE
HELIXER_GFF is the window-coordinate GFF3 cache: reused when it exists, written otherwise.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import panel_helixer_metric as pm  # noqa: E402
import yaml  # noqa: E402
from Bio import SeqIO  # noqa: E402
from MATPredict.detect import tracks as T  # noqa: E402
from MATPredict.detect.helixer_rule import REFINED, refine_gene  # noqa: E402
from MATPredict.detect.helixer_windows import merged_windows, write_window_fasta  # noqa: E402
from MATPredict.detect.reference_fasta import build_reference_fasta  # noqa: E402


def locus_spans(doc: dict):
    return [(d["contig"], int(d["start"]), int(d["end"])) for d in doc.get("detected") or []
            if d.get("contig") and d.get("start") is not None and d.get("end") is not None]


def query_index(db_root: Path) -> dict[tuple[str, str], str]:
    """(record_id, gene name) -> curated protein, from every record in the database."""
    with tempfile.TemporaryDirectory() as tmp:
        fa = build_reference_fasta(db_root, Path(tmp) / "ref.faa")
        out = {}
        for rec in SeqIO.parse(str(fa), "fasta"):
            record_id, _, name = rec.id.split("|", 2)
            out[(record_id, name)] = str(rec.seq)
        return out


def protein(exons, strand, seq: str) -> str:
    from Bio.Seq import Seq
    dna = T.coding_sequence(exons, strand, seq)
    return str(Seq(dna[: len(dna) // 3 * 3]).translate())


def helixer_track(genome: Path, spans, gff: Path, exe: str, cpus: int, seqs, lens) -> T.Track:
    if not gff.exists():
        wins = merged_windows(spans, 5000, lens)
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / "w.fa"
            write_window_fasta(wins, seqs, fa)
            subprocess.run([exe, "-f", str(fa), "-o", str(gff), "-l", "fungi", "-c", str(cpus)],
                           check=True, capture_output=True)
    return T.from_helixer_gff3(gff, seqs, window_offsets=True)


def check_one(name, genome, detect_dir, ncbi_gff, hx_gff, queries, exe, cpus, totals, extent="gene") -> None:
    doc = yaml.safe_load((Path(detect_dir) / "detection_report.yaml").read_text())
    spans = locus_spans(doc)
    seqs, lens = {}, {}
    wanted = {c for c, _, _ in spans}
    for r in SeqIO.parse(str(genome), "fasta"):
        lens[r.id] = len(r.seq)
        if r.id in wanted:
            seqs[r.id] = str(r.seq).upper()
    track = helixer_track(Path(genome), spans, Path(hx_gff), exe, cpus, seqs, lens)
    ncbi = pm.read_ncbi_genes(Path(ncbi_gff), extent)
    print(f"== {name}")
    for call in doc.get("detected") or []:
        for g in call.get("gene_evidence") or []:
            if not (g.get("role") or "").startswith("flanking") or g.get("lost_to") or g.get("start") is None:
                continue
            o = (g["contig"], g["start"], g["end"], g["strand"], g["gene"])
            n = pm.best(o, ncbi)
            if n is None:
                continue
            exons = tuple(sorted((x["start"], x["end"]) for x in g["exons"])) if g.get("exons") else \
                ((g["start"], g["end"]),)
            gene = SimpleNamespace(gene_name=g["gene"], role=g["role"], contig=g["contig"], start=g["start"],
                                   end=g["end"], strand=g["strand"])
            q = queries.get((g.get("reference_record"), g["gene"]))
            unver = refine_gene(gene, track)
            ver = refine_gene(gene, track, query=q, homology_protein=protein(exons, g["strand"], seqs[g["contig"]]),
                              helixer_protein=lambda h: protein(h.exons, h.strand, seqs[h.contig])) \
                if q else unver
            base = pm.within(n, o)
            u_ok = pm.within(n, (g["contig"], unver.start, unver.end, g["strand"]))
            v_ok = pm.within(n, (g["contig"], ver.start, ver.end, g["strand"]))
            for key, ok, r in (("unverified", u_ok, unver), ("verified", v_ok, ver)):
                totals[key]["n"] += 1
                totals[key]["bounded"] += ok
                totals[key]["changed"] += r.action == REFINED
                totals[key]["gained"] += ok and not base
                totals[key]["lost"] += base and not ok
            totals["baseline"]["n"] += 1
            totals["baseline"]["bounded"] += base
            gain = "n/a" if ver.score_gain is None else f"{ver.score_gain:+.0f}"
            print(f"   {g['gene']:6s} {g['strand']} NCBI {n[1]}-{n[2]} | homology {'ok ' if base else 'off'} "
                  f"| unverified {unver.action:7s} {'ok ' if u_ok else 'off'} | query gain {gain:>5s} "
                  f"-> {ver.action:7s} {'ok ' if v_ok else 'off'}  {ver.reason}"
                  + ("" if q else "  [no query found]"))


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db-root", required=True)
    p.add_argument("--helixerlite", required=True)
    p.add_argument("--cpus", type=int, default=4)
    p.add_argument("--ncbi-extent", choices=["gene", "cds"], default="gene",
                   help="compare with the NCBI gene feature (default) or the span of its CDS features")
    p.add_argument("--run", nargs=5, action="append", required=True,
                   metavar=("NAME", "GENOME_FNA", "DETECT_DIR", "NCBI_GFF", "HELIXER_GFF"))
    a = p.parse_args(argv)
    queries = query_index(Path(a.db_root))
    totals = {k: {"n": 0, "bounded": 0, "changed": 0, "gained": 0, "lost": 0}
              for k in ("baseline", "unverified", "verified")}
    for name, genome, detect_dir, ncbi, hx in a.run:
        check_one(name, genome, detect_dir, ncbi, hx, queries, a.helixerlite, a.cpus, totals, a.ncbi_extent)
    print("\nflank genes matching NCBI, ends within 30 bp:")
    for k, t in totals.items():
        extra = "" if k == "baseline" else f"  changed {t['changed']}, gained {t['gained']}, lost {t['lost']}"
        print(f"  {k:11s} {t['bounded']}/{t['n']}{extra}")


if __name__ == "__main__":
    main()
