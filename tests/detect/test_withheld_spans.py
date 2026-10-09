"""Withheld loci report `core_span` and `supported_span` too (so the span question can be asked of the 164-style changes).

Seam: run_pipeline (injected search) and the report writer. A locus withheld at the fraction floor is the fixture:
2 of 5 expected genes (0.4 < 0.5), so it is not called but is listed in `suppressed_loci`.
"""
from __future__ import annotations

import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import BELOW_FRACTION_FLOOR, run_pipeline
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _model
from tests.detect.test_below_fraction_floor_reported import _order

KEY = FamilyKey("P", "MAT")


def _hit(gene, role, start, end, bits, identity=100.0):
    return SearchHit(KEY, gene, role, "c1", start, end, "+", identity, "rec1", "tblastn_genome", bitscore=bits, evalue=1e-20 if bits > 39 else 4.0)


def _outcome(tmp_path):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(_order(4))
    rec = tmp_path / "P" / "Fam" / "rec1"
    rec.mkdir(parents=True)
    rec.joinpath("metadata.yaml").write_text("record_id: rec1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n")
    hits = [_hit("sexP", "core_MAT", 1000, 1600, 120.0), _hit("f1", "flanking_variable", 3000, 4000, 100.0, 99.0),
            _hit("f1", "flanking_variable", 20000, 20300, 90.0, 40.0),     # strong extra copy, not a reported gene
            _hit("f1", "flanking_variable", 28000, 28100, 25.0, 40.0)]     # weak extra copy
    spans = {"sexP": (1000, 1600)}
    return run_pipeline(
        genome_fasta=tmp_path / "g.fa", proteome_fasta=None, taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "r.faa",
        search_localize=lambda *a, **k: hits,
        polish_with_exonerate=lambda *, gene_name, **k: (_model(gene_name, "c1", *spans[gene_name], family_key=KEY, role="core_MAT",
                                                              method="exonerate_refine") if gene_name in spans else None),
        polish_with_miniprot=lambda *, gene_name, **k: None, relaxed_second_pass=False)


def test_a_withheld_locus_reports_its_core_and_supported_spans(tmp_path):
    outcome = _outcome(tmp_path)
    out = tmp_path / "r.yaml"
    write_detection_report(outcome, out)
    (row,) = [s for s in yaml.safe_load(out.read_text())["suppressed_loci"] if s["withheld_reason"] == BELOW_FRACTION_FLOOR]
    assert (row["start"], row["end"]) == (1000, 28100)                                      # cluster span, as before
    assert row["core_span"] == {"start": 1000, "end": 4000, "beyond_core_bp": 24100}
    assert row["supported_span"] == {"start": 1000, "end": 20300, "min_bitscore": 39.0, "beyond_supported_bp": 7800}
