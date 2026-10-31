"""A GenBank genome with sequence but no CDS translations is searched as a genome, not refused (curator, 2026-10-30)."""
import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from MATPredict.detect.cli import _use_genbank_genome
from MATPredict.detect.genbank_input import prepare_genbank_input

CDS = "ATGGCTGCTAAAGGTCCAGATTTCTTGGAAGGAGCTGCTAGAAGCCCATATGCAGAAGCTTAA"
CONTIG = "N" * 100 + CDS + "N" * 100


def _write(path: Path, translation: bool) -> Path:
    rec = SeqRecord(Seq(CONTIG), id="c1.1", name="c1", description="t")
    rec.annotations["molecule_type"] = "DNA"
    q = {"locus_tag": ["X_1"]}
    if translation:
        q["translation"] = [str(Seq(CDS).translate(to_stop=True))]
    rec.features.append(SeqFeature(FeatureLocation(100, 163, strand=1), type="CDS", qualifiers=q))
    SeqIO.write(rec, path, "genbank")
    return path


def test_sequence_only_genbank_writes_the_genome_and_no_proteome(tmp_path):
    gb = prepare_genbank_input(_write(tmp_path / "a.gbk", translation=False), tmp_path / "w")
    assert gb.converted.proteome_fasta is None and gb.converted.proteins == 0 and gb.converted.contigs == 1
    assert gb.converted.genome_fasta.read_text().startswith(">c1.1")


def test_detect_searches_the_genome_of_a_sequence_only_genbank(tmp_path):
    args = argparse.Namespace(genome=str(_write(tmp_path / "a.gbk", translation=False)), proteins=None, annotation=None)
    _use_genbank_genome(args, tmp_path / "out")
    assert args.genome.endswith("a.fna") and args.proteins is None and args.annotation is None


def test_annotated_genbank_still_takes_the_proteome_route(tmp_path):
    src = _write(tmp_path / "b.gbk", translation=True)
    args = argparse.Namespace(genome=str(src), proteins=None, annotation=None)
    _use_genbank_genome(args, tmp_path / "out")
    assert args.proteins.endswith("b.faa") and args.annotation == str(src)
