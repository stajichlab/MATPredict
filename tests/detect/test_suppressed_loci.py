"""A withheld locus keeps its coordinates, and its not-detected reason is true.

The modelled-gene bar withholds loci, and before this the report kept only a
count. That made bar losses unauditable: on Pezizales, 81 of 172 genomes were
withheld by the bar and no report said WHERE the withheld cluster sat, so it
could not be checked against a known locus. The holdout scorer needs the same
coordinates to tell a bar loss from a genuine miss.
"""
import yaml

from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.report import write_detection_report

from tests.detect.test_modelled_gene_bar import _annotated, _raw, _run
from tests.detect.test_pipeline import FAMILY, _no_polish, _write_order, _write_record


def test_a_withheld_locus_keeps_its_coordinates(tmp_path):
    outcome = _run(tmp_path, [_annotated("mfa1", 100, 200), _raw("pra1", 300, 400)])
    assert outcome.results == []
    assert len(outcome.suppressed_loci) == 1
    locus = outcome.suppressed_loci[0]
    assert (locus.contig, locus.start, locus.end) == ("c1", 100, 400)
    assert locus.polished_genes == 1


def test_the_report_writes_the_withheld_loci(tmp_path):
    outcome = _run(tmp_path, [_annotated("mfa1", 100, 200), _raw("pra1", 300, 400)])
    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["suppressed_unpolished"] == 1
    [entry] = doc["suppressed_loci"]
    assert entry["contig"] == "c1"
    assert (entry["start"], entry["end"]) == (100, 400)
    assert entry["polished_genes"] == 1
    assert entry["family"] == "P:aLocus"
    assert set(entry["genes_found"]) == {"mfa1", "pra1"}


def test_the_reason_states_the_real_modelled_count(tmp_path):
    """One gene WAS modelled; 'never modelled' would be false."""
    outcome = _run(tmp_path, [_annotated("mfa1", 100, 200), _raw("pra1", 300, 400)])
    [reason] = [n.reason for n in outcome.not_detected if n.family_key == FAMILY.key]
    assert "never modelled" not in reason
    assert "1 modelled gene" in reason


def test_a_zero_model_reason_still_says_none_could_be_modelled(tmp_path):
    outcome = _run(tmp_path, [_raw("mfa1", 100, 200), _raw("pra1", 300, 400)])
    [reason] = [n.reason for n in outcome.not_detected if n.family_key == FAMILY.key]
    assert "none" in reason


def test_a_family_with_no_reference_proteins_says_so(tmp_path):
    """A holdout can empty a family. 'No hits' would blame the genome."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    (tmp_path / "reference.faa").write_text("")
    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=lambda *a, **k: [],
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    [reason] = [n.reason for n in outcome.not_detected if n.family_key == FAMILY.key]
    assert "no reference protein" in reason


def test_a_failed_genetic_code_lookup_is_recorded(tmp_path, monkeypatch):
    """Serinales translate CTG as serine (table 12). A failed lookup falls back
    to table 1; that must be visible, since it changes every translation."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def boom(_taxid):
        raise RuntimeError("429 Too Many Requests")

    from MATPredict.detect import pipeline
    from MATPredict.detect.family_registry import RoutingDecision
    monkeypatch.setattr(pipeline, "route", lambda taxid, fams, **k: RoutingDecision(fams, "lineage"))
    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=5476,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=lambda *a, **k: [],
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        genetic_code_resolver=boom,
    )
    assert outcome.genetic_code == 1
    assert "429" in outcome.genetic_code_error
    path = tmp_path / "r.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["genetic_code"] == 1 and "429" in doc["genetic_code_error"]
