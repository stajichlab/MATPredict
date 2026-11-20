"""Fixture test of scripts/model_vs_annotation.py (#290, #349): intron matching against a GFF3 and the end-to-end verdict."""
from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

import yaml

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("mva", SCRIPTS / "model_vs_annotation.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

GENOME = "CCCC" + "ATGGCTGCTAAA" + "GT" + "A" * 40 + "AG" + "GCTAAAGGTTAA" + "CCCC"
GOOD = [(5, 16), (61, 72)]
SHIFTED = [(5, 18), (63, 72)]
GFF = ("c1\tx\tmRNA\t5\t72\t.\t+\t.\tID=m1\n" "c1\tx\tCDS\t5\t16\t.\t+\t0\tID=m1-CDS1;Parent=m1\n" "c1\tx\tCDS\t61\t72\t.\t+\t0\tID=m1-CDS2;Parent=m1\n")


def write_gff(path: Path):
    with gzip.open(path, "wt") as fh:
        fh.write("##gff-version 3\n" + GFF)


def test_load_and_match(tmp_path):
    write_gff(tmp_path / "a.gff.gz")
    m = mod.load_mrnas(tmp_path / "a.gff.gz")["c1"][0]
    assert m["introns"] == [(17, 60)] and m["strand"] == "+"
    assert mod.matches(GOOD, m) and not mod.matches(SHIFTED, m)
    assert mod.matches([(5, 16)], m)                       # a model that stops before the intron has no intron in its span
    assert mod.best_mrna({"c1": [m]}, "c1", "+", GOOD) is m and mod.best_mrna({"c1": [m]}, "c1", "-", GOOD) is None
    assert mod.verdict(True, False) == "primary_matches" and mod.verdict(False, False) == "neither_matches"


def entry(primary, alternate):
    ex = lambda es: [{"start": s, "end": e} for s, e in es]
    return {"gene": "G1", "role": "core_MAT", "status": "polished_disagree", "contig": "c1", "start": 5, "end": 72, "strand": "+", "method": "exonerate_refine",
            "identity": 99.0, "exons": ex(primary), "alternate_model": {"contig": "c1", "start": 5, "end": 72, "strand": "+", "method": "miniprot_refine", "identity": 90.0, "exons": ex(alternate)}}


def doc(e):
    return {"genetic_code": 1, "detected": [{"family": "Ascomycota:MAT", "confidence": "high", "start": 5, "end": 72, "gene_evidence": [e]}]}


def test_end_to_end_verdicts(tmp_path):
    for s in ("A", "B", "N"):
        with gzip.open(tmp_path / f"{s}_g.fa.gz", "wt") as fh:
            fh.write(">c1\n" + GENOME + "\n")
    write_gff(tmp_path / "A.gff.gz")
    write_gff(tmp_path / "B.gff.gz")
    docs = [("A", doc(entry(SHIFTED, GOOD))), ("B", doc(entry(GOOD, SHIFTED))), ("N", doc(entry(GOOD, SHIFTED)))]   # N has no GFF
    rows = mod.run(docs, str(tmp_path / "{sample}_g.fa.gz"), str(tmp_path / "{sample}.gff.gz"), tmp_path / "db", tmp_path / "o" / "t")
    by = {r["sample"]: r for r in rows}
    assert by["A"]["annotation"] == "alternate_matches" and by["A"]["class"] == "primary_issue_alt_clean" and by["A"]["trigger_noncanonical"] is True
    assert by["B"]["annotation"] == "primary_matches" and by["B"]["class"] == "primary_clean_alt_issue"
    assert by["N"]["annotation"] == "no_gff"
    # an annotated gene on the other strand overlaps nothing: no_overlap
    write_gff(tmp_path / "S.gff.gz")
    with gzip.open(tmp_path / "S_g.fa.gz", "wt") as fh:
        fh.write(">c1\n" + GENOME + "\n")
    e = entry(GOOD, SHIFTED)
    e["strand"] = e["alternate_model"]["strand"] = "-"
    assert mod.run([("S", doc(e))], str(tmp_path / "{sample}_g.fa.gz"), str(tmp_path / "{sample}.gff.gz"), tmp_path / "db", tmp_path / "o" / "s")[0]["annotation"] in ("no_overlap", "not_checked")
    assert "all\tprimary_issue_alt_clean\talternate_matches\t1" in (tmp_path / "o" / "t.summary.tsv").read_text()
