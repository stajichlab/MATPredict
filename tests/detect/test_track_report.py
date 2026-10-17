"""Gene-model tracks in the report, and GenBank as a detect input."""
import sys
from argparse import Namespace
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from MATPredict.detect import cli
from MATPredict.detect.genbank_input import is_genbank
from MATPredict.detect.track_report import build_tracks_doc, load_tracks
from MATPredict.report.html import render_genome_report

CDS = "ATGGCTGCTAAAGGTCCAGATTTCTTGGAAGGAGCTGCTAGAAGCCCATATGCAGAAGCTTAA"
CONTIG = "N" * 100 + CDS + "N" * 100          # CDS at 101-160


def _write_genbank(path: Path, with_sequence: bool = True, end: int = 163, translation: bool = True) -> Path:
    rec = SeqRecord(Seq(CONTIG) if with_sequence else Seq(None, length=len(CONTIG)), id="c1", name="c1",
                    description="t")
    rec.annotations["molecule_type"] = "DNA"
    q = {"gene": ["sexM"], "locus_tag": ["X_1"]}
    if translation:
        q["translation"] = [str(Seq(CDS).translate(to_stop=True))]
    rec.features.append(SeqFeature(FeatureLocation(100, end, strand=1), type="CDS", qualifiers=q))
    SeqIO.write(rec, path, "genbank")
    return path


def _genome(tmp_path: Path) -> Path:
    fa = tmp_path / "g.fna"
    fa.write_text(f">c1\n{CONTIG}\n")
    return fa


def _report_doc() -> dict:
    return {"detected": [{
        "contig": "c1", "start": 101, "end": 163,
        "gene_evidence": [{"gene": "sexM", "role": "core_MAT", "contig": "c1", "start": 101, "end": 163,
                           "strand": "+", "exons": [{"start": 101, "end": 163}]}]}]}


def test_is_genbank_by_suffix_and_by_content(tmp_path: Path):
    gbk = _write_genbank(tmp_path / "a.txt")
    assert is_genbank(gbk) and is_genbank(tmp_path / "x.gbff")
    assert not is_genbank(_genome(tmp_path))


def test_genbank_track_agrees_with_the_homology_model(tmp_path: Path):
    doc = build_tracks_doc(_report_doc(), _genome(tmp_path), genbank=_write_genbank(tmp_path / "a.gbk"))
    (row,) = doc["summaries"]["genbank"]
    assert row["gene"] == "sexM" and row["model"]["other_gene"] == "sexM"
    assert row["model"]["ends_agree"] and row["model"]["other_complete_orf"] is True
    assert doc["sources"]["evidence"]["n_models"] == 1 and doc["sources"]["genbank"]["n_models"] == 1


def test_genbank_model_with_longer_end_is_reported_as_extension(tmp_path: Path):
    gbk = _write_genbank(tmp_path / "a.gbk", end=190)
    (row,) = build_tracks_doc(_report_doc(), _genome(tmp_path), genbank=gbk)["summaries"]["genbank"]
    assert (row["model"]["extends_5p"], row["model"]["extends_3p"]) == (0, 27)
    assert row["model"]["ends_agree"]          # within the 30 bp tolerance


def test_nothing_to_compare_gives_no_tracks(tmp_path: Path):
    assert build_tracks_doc(_report_doc(), _genome(tmp_path)) is None
    assert build_tracks_doc({"detected": []}, _genome(tmp_path), genbank=_write_genbank(tmp_path / "a.gbk")) is None


FAKE_HELIXER = '''
import sys
a = sys.argv[1:]
fasta, out = a[a.index("-f") + 1], a[a.index("-o") + 1]
name = open(fasta).readline()[1:].strip()
open(out, "w").write(
    "##gff-version 3\\n"
    f"{name}\\tHelixer\\tgene\\t1\\t60\\t.\\t+\\t.\\tID=g1\\n"
    f"{name}\\tHelixer\\tmRNA\\t1\\t60\\t.\\t+\\t.\\tID=g1.t;Parent=g1\\n"
    f"{name}\\tHelixer\\tCDS\\t1\\t60\\t.\\t+\\t0\\tID=g1.c;Parent=g1.t\\n")
'''


def test_helixer_track_from_a_command_on_the_locus_window(tmp_path: Path):
    fake = tmp_path / "fake_helixer.py"
    fake.write_text(FAKE_HELIXER)
    # window is the whole 260 bp contig (clipped), so window-relative 1-60 is genome 1-60: no overlap with 101-160
    doc = build_tracks_doc(_report_doc(), _genome(tmp_path), helixer_cmd=[sys.executable, str(fake)], flank=5000)
    assert doc["sources"]["helixer"]["n_models"] == 1
    (row,) = doc["summaries"]["helixer"]
    assert row["model"] is None                 # reported as 'no model', not hidden


def test_html_section_shows_tracks_and_is_absent_without_them(tmp_path: Path):
    tracks = build_tracks_doc(_report_doc(), _genome(tmp_path), genbank=_write_genbank(tmp_path / "a.gbk"))
    html = render_genome_report({"detected": [], "run": {}}, tracks=tracks, generated="t")
    assert "Gene-model tracks" in html and "GenBank annotation" in html and "1 of 1 homology gene" in html
    assert "Gene-model tracks" not in render_genome_report({"detected": [], "run": {}}, generated="t")


def test_cli_genbank_genome_becomes_fasta_proteome_and_annotation(tmp_path: Path):
    gbk = _write_genbank(tmp_path / "a.gbk")
    args = Namespace(genome=str(gbk), proteins=None, annotation=None)
    cli._use_genbank_genome(args, tmp_path / "out")
    assert args.genome.endswith("a.fna") and args.proteins.endswith("a.faa") and args.annotation == str(gbk)
    assert ">X_1 c1:101-163:+" in Path(args.proteins).read_text()


def test_cli_keeps_a_given_proteome_and_ignores_fasta(tmp_path: Path):
    args = Namespace(genome=str(_genome(tmp_path)), proteins="mine.faa", annotation=None)
    cli._use_genbank_genome(args, tmp_path / "out")
    assert args.proteins == "mine.faa" and args.annotation is None
    gbk = _write_genbank(tmp_path / "a.gbk")
    args = Namespace(genome=str(gbk), proteins="mine.faa", annotation=None)
    cli._use_genbank_genome(args, tmp_path / "out")
    assert args.proteins == "mine.faa" and args.genome.endswith("a.fna")


def test_genbank_without_sequence_is_refused_with_advice(tmp_path: Path):
    import pytest
    gbk = _write_genbank(tmp_path / "a.gbk", with_sequence=False)
    with pytest.raises(ValueError, match="--annotation"):
        cli._use_genbank_genome(Namespace(genome=str(gbk), proteins=None, annotation=None), tmp_path / "out")


def test_write_tracks_step_writes_yaml_and_report_reads_it(tmp_path: Path):
    out = tmp_path / "run"
    out.mkdir()
    (out / "detection_report.yaml").write_text(yaml.safe_dump(_report_doc()))
    gbk = _write_genbank(tmp_path / "a.gbk")
    cli._write_tracks(Namespace(genome=str(_genome(tmp_path)), annotation=str(gbk), helixerlite=None), out)
    assert load_tracks(out)["sources"]["genbank"]["n_models"] == 1
    cli._write_tracks(Namespace(genome=str(_genome(tmp_path)), annotation=None, helixerlite=None),
                      tmp_path / "other")           # no sources: nothing written, no error
    assert not (tmp_path / "other").exists()
