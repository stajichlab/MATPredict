"""Fixture test of scripts/cap_annotation_review.py (#341)."""
from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("car", Path(__file__).resolve().parents[1] / "scripts" / "cap_annotation_review.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

GENES = {"c1": [(100, 400), (1000, 1100)]}


def test_verdicts():
    assert mod.verdict(GENES, "c1", 150, 350) == "annotated"
    assert mod.verdict(GENES, "c1", 300, 900) == "partial"       # 101 of 601 bases
    assert mod.verdict(GENES, "c1", 500, 900) == "novel"
    assert mod.verdict(GENES, "c2", 1, 10) == "novel"


def test_folder_candidates():
    assert mod.folder_candidates("Aaosphaeria arxii CBS 175.79", "CBS 175.79") == ["Aaosphaeria_arxii_CBS_175.79"]
    assert mod.folder_candidates("Corynespora cassiicola", "IA") == ["Corynespora_cassiicola_IA", "Corynespora_cassiicola"]


def ev(gene, status, bits, start, end):
    return {"gene": gene, "status": status, "bitscore": bits, "contig": "c1", "start": start, "end": end}


def test_review_ranks_only_the_called_idiomorph_and_classifies():
    locus = {"idiomorph": "MAT1-1", "gene_evidence": [ev("MAT1-1-1", "polished_single", 41.0, 150, 350), ev("MAT1-2-1", "polished_agree", 300.0, 500, 900)],
             "genes_found": ["MAT1-1-1", "MAT1-2-1"]}
    r = mod.review_locus(locus, GENES)
    assert r["strongest_own_gene"] == "MAT1-1-1" and r["review"] == "annotated_gene_present"
    assert "MAT1-2-1:polished_agree:300.0:novel" in r["other_idiomorph_genes"]
    assert mod.review_locus(locus, None)["review"] == "no_annotation"
    locus["gene_evidence"][0]["start"], locus["gene_evidence"][0]["end"] = 500, 900
    assert mod.review_locus(locus, GENES)["review"] == "no_annotated_gene"


def test_run_end_to_end(tmp_path):
    ann = tmp_path / "ann" / "Sp_a_S1" / "predict_results"
    ann.mkdir(parents=True)
    with gzip.open(ann / "x.gff3.gz", "wt") as fh:
        fh.write("c1\tfunannotate\tgene\t100\t400\t.\t+\t.\tID=g1\n")
    (tmp_path / "samples.csv").write_text("ASMID,SPECIES_IN,STRAIN\nG1,Sp a,S1\nG2,Sp b,S2\n")
    (tmp_path / "changes.tsv").write_text("genome\tcontig\tstart\tend\tnew_rule\nG1\tc1\t100\t400\tcapped\nG2\tc1\t100\t400\tcapped\n")
    mk = lambda: {"detected": [{"family": "Ascomycota:MAT", "contig": "c1", "start": 100, "end": 400, "idiomorph": "MAT1-1", "confidence": "high",
                                "gene_evidence": [ev("MAT1-1-1", "unpolished", 41.0, 150, 350)], "genes_found": ["MAT1-1-1"]}]}
    counts = mod.run([("G1", mk()), ("G2", mk())], tmp_path / "changes.tsv", tmp_path / "samples.csv", tmp_path / "ann", tmp_path / "out")
    assert counts == {"annotated_gene_present": 1, "no_annotation": 1}
    assert (tmp_path / "out" / "review.tsv").read_text().count("\n") == 3
