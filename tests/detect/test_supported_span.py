"""`supported_span`: the reported span rebuilt from modelled genes plus hits strong enough to trust (report only).

The cluster span is chained from every hit at BLAST's default e-value, so one weak hit can stretch a locus. The supported
span keeps the locus's own modelled genes and its own hits at or above a bitscore floor, and drops the rest. The cluster
span and every call stay as they are. Seam: run_pipeline (injected search) and the report writers.
"""
from __future__ import annotations

import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.report import write_detection_gff3, write_detection_report
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _model

KEY = FamilyKey("P", "MAT")
ORDER = (
    "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n    idiomorph_values: [Plus, Minus]\n"
    "    taxonomic_scope: [1]\n    genes:\n      - {name: f1, role: flanking_variable}\n"
    "      - {name: f2, role: flanking_variable}\n      - {name: sexP, role: core_MAT, present_in_idiomorphs: [Plus]}\n"
)


def _hit(gene, role, start, end, bits, ev, identity=60.0):
    return SearchHit(KEY, gene, role, "c1", start, end, "+", identity, "rec1", "tblastn_genome", bitscore=bits, evalue=ev)


def _run(tmp_path, extra_hits, **kw):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(ORDER)
    rec = tmp_path / "P" / "Fam" / "rec1"
    rec.mkdir(parents=True)
    rec.joinpath("metadata.yaml").write_text("record_id: rec1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n")
    hits = [_hit("sexP", "core_MAT", 1000, 1600, 120.0, 1e-30), _hit("f1", "flanking_variable", 3000, 4000, 100.0, 1e-25)] + extra_hits
    spans = {"sexP": (1000, 1600), "f1": (3000, 4000)}

    def model(gene, method):
        return (_model(gene, "c1", *spans[gene], family_key=KEY, role="core_MAT" if gene == "sexP" else "flanking_variable",
                       method=method) if gene in spans else None)

    return run_pipeline(
        genome_fasta=tmp_path / "g.fa", proteome_fasta=None, taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "r.faa",
        search_localize=lambda *a, **k: hits,
        polish_with_exonerate=lambda *, gene_name, **k: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **k: model(gene_name, "miniprot_refine"),
        relaxed_second_pass=False, **kw)


# Extra copies of the already-reported gene f1 (lower identity than its best hit, so they never become the gene's evidence):
# gene_evidence holds the best hit of every gene the locus reports, and those always count.
STRONG = _hit("f1", "flanking_variable", 20000, 20300, 90.0, 1e-18, identity=40.0)   # own family, strong, not a reported gene
WEAK = _hit("f1", "flanking_variable", 28000, 28100, 25.0, 4.0, identity=40.0)       # the T48-F kind of hit


def test_a_weak_hit_stretches_the_cluster_span_but_not_the_supported_span(tmp_path):
    outcome = _run(tmp_path, [STRONG, WEAK])
    (r,) = outcome.results
    assert (r.start, r.end) == (1000, 28100)                       # cluster span, unchanged behaviour
    assert r.supported_span == {"start": 1000, "end": 20300, "min_bitscore": 39.0, "beyond_supported_bp": 7800}


def test_the_floor_is_a_parameter(tmp_path):
    low, high = tmp_path / "low", tmp_path / "high"
    low.mkdir(); high.mkdir()
    (r,) = _run(low, [STRONG, WEAK], supported_min_bitscore=20.0).results
    assert r.supported_span["end"] == 28100 and r.supported_span["beyond_supported_bp"] == 0
    (r,) = _run(high, [STRONG, WEAK], supported_min_bitscore=95.0).results
    assert r.supported_span["end"] == 4000                          # 90 bits no longer passes; the modelled genes always count


def test_reported_genes_always_count_whatever_their_bitscore(tmp_path):
    hits_low = [STRONG]
    (r,) = _run(tmp_path, hits_low, supported_min_bitscore=1000.0).results
    assert (r.supported_span["start"], r.supported_span["end"]) == (1000, 4000)


def test_the_report_and_gff3_carry_the_supported_span(tmp_path):
    outcome = _run(tmp_path, [STRONG, WEAK])
    y = tmp_path / "r.yaml"; g = tmp_path / "r.gff3"
    write_detection_report(outcome, y); write_detection_gff3(outcome, g)
    locus = yaml.safe_load(y.read_text())["detected"][0]
    assert locus["supported_span"] == {"start": 1000, "end": 20300, "min_bitscore": 39.0, "beyond_supported_bp": 7800}
    assert (locus["start"], locus["end"]) == (1000, 28100) and locus["core_span"]["end"] == 4000
    line = next(l for l in g.read_text().splitlines() if "\tMAT_locus\t" in l)
    attrs = dict(kv.split("=", 1) for kv in line.split("\t")[8].split(";") if "=" in kv)
    assert (attrs["supported_start"], attrs["supported_end"], attrs["beyond_supported_bp"]) == ("1000", "20300", "7800")
    assert (line.split("\t")[3], line.split("\t")[4]) == ("1000", "28100")


def test_a_hit_without_a_bitscore_is_not_trusted(tmp_path):
    unknown = SearchHit(KEY, "f1", "flanking_variable", "c1", 20000, 20300, "+", 40.0, "rec1", "tblastn_genome")
    (r,) = _run(tmp_path, [unknown]).results
    assert r.supported_span["end"] == 4000


def test_a_strong_hit_of_another_family_stretches_the_cluster_but_is_not_support(tmp_path):
    other = FamilyKey("P", "OTHER")
    foreign = SearchHit(other, "g1", "core_MAT", "c1", 20000, 20300, "+", 60.0, "rec2", "tblastn_genome", bitscore=100.0, evalue=1e-25)
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(ORDER + (
        "  - locus_name: OTHER\n    vocabulary_type: enum\n    idiomorph_values: [A, B]\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: g1, role: core_MAT}\n"))
    for rid, loc in (("rec1", "MAT"), ("rec2", "OTHER")):
        rec = tmp_path / "P" / "Fam" / rid; rec.mkdir(parents=True)
        rec.joinpath("metadata.yaml").write_text(f"record_id: {rid}\nmating_type: {{locus_name: {loc}, idiomorphs: [{'Plus' if loc == 'MAT' else 'A'}]}}\n")
    spans = {"sexP": (1000, 1600), "f1": (3000, 4000)}
    hits = [_hit("sexP", "core_MAT", 1000, 1600, 120.0, 1e-30), _hit("f1", "flanking_variable", 3000, 4000, 100.0, 1e-25), foreign]
    outcome = run_pipeline(
        genome_fasta=tmp_path / "g.fa", proteome_fasta=None, taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "r.faa",
        search_localize=lambda *a, **k: hits,
        polish_with_exonerate=lambda *, gene_name, **k: _model(gene_name, "c1", *spans[gene_name], family_key=KEY,
            role="core_MAT" if gene_name == "sexP" else "flanking_variable", method="exonerate_refine") if gene_name in spans else None,
        polish_with_miniprot=lambda *, gene_name, **k: None, relaxed_second_pass=False)
    (r,) = [x for x in outcome.results if x.family_key == KEY]
    assert (r.start, r.end) == (1000, 20300)            # the cluster span includes the other family's hit (the defect)
    assert r.supported_span["end"] == 4000 and r.supported_span["beyond_supported_bp"] == 16300


def test_cli_option_defaults_to_the_flank_carried_floor_and_parses():
    from MATPredict.__main__ import build_parser
    from MATPredict.detect.family_registry import DEFAULT_FLANK_CARRIED_MIN_BITSCORE
    p = build_parser()
    base = ["detect", "--genome", "g.fa", "--out-dir", "o"]
    assert p.parse_args(base).supported_min_bitscore == DEFAULT_FLANK_CARRIED_MIN_BITSCORE == 39.0
    assert p.parse_args(base + ["--supported-min-bitscore", "50"]).supported_min_bitscore == 50.0
