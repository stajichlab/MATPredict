"""The genome FASTA is indexed once per genome, not once per polish call.

`_extract_window` is called twice for every gene polished (exonerate and
miniprot each get the same window), and once more for every additional gene
and every additional admitted family in a cluster. Each call used to run
`SeqIO.index` over the WHOLE genome. Measured on a synthetic 46 MB assembly
(300 contigs, the shape of an ordinary Agaricomycete): 565 ms per call, so a
cluster with 10 genes admitted to 3 families spent ~34 s re-reading the same
file before any alignment ran.
"""
from pathlib import Path

import pytest
from Bio import SeqIO

from MATPredict.detect import search as search_module
from MATPredict.detect.search import _extract_window


def _genome(tmp_path: Path) -> Path:
    path = tmp_path / "genome.fa"
    path.write_text(">c1\n" + "ACGT" * 250 + "\n>c2\n" + "TTTT" * 250 + "\n")
    return path


@pytest.fixture(autouse=True)
def _clear_cache():
    search_module.clear_genome_index_cache()
    yield
    search_module.clear_genome_index_cache()


def _counting_index(monkeypatch) -> list[int]:
    calls: list[int] = []
    real = SeqIO.index

    def counted(path, fmt, *a, **kw):
        calls.append(1)
        return real(path, fmt, *a, **kw)

    monkeypatch.setattr(search_module.SeqIO, "index", counted)
    return calls


def test_repeated_windows_on_one_genome_index_it_once(tmp_path, monkeypatch):
    genome = _genome(tmp_path)
    calls = _counting_index(monkeypatch)
    for start in (1, 5, 9, 13):
        _extract_window(genome, ("c1", start, start + 20), tmp_path)
    assert len(calls) == 1


def test_the_extracted_sequence_is_still_correct(tmp_path, monkeypatch):
    """Caching must not change what comes out. Same slice, twice, and it must
    equal the real subsequence both times."""
    genome = _genome(tmp_path)
    expected = str(next(SeqIO.parse(genome, "fasta")).seq)[9:29]
    for _ in range(2):
        out = _extract_window(genome, ("c1", 10, 29), tmp_path)
        record = next(SeqIO.parse(out, "fasta"))
        assert str(record.seq) == expected
        assert record.id == "c1"


def test_a_second_genome_is_indexed_separately(tmp_path, monkeypatch):
    first = _genome(tmp_path)
    second = tmp_path / "other.fa"
    second.write_text(">c1\n" + "GGGG" * 250 + "\n")
    calls = _counting_index(monkeypatch)
    _extract_window(first, ("c1", 1, 20), tmp_path)
    _extract_window(second, ("c1", 1, 20), tmp_path)
    _extract_window(first, ("c1", 30, 50), tmp_path)
    assert len(calls) == 2


def test_a_rewritten_genome_is_reindexed(tmp_path, monkeypatch):
    """The cache is keyed on the file's size and mtime, so a genome replaced
    on disk mid-session is never served from a stale index."""
    genome = _genome(tmp_path)
    calls = _counting_index(monkeypatch)
    _extract_window(genome, ("c1", 1, 20), tmp_path)
    genome.write_text(">c1\n" + "GGGGCCCC" * 200 + "\n")
    out = _extract_window(genome, ("c1", 1, 8), tmp_path)
    assert len(calls) == 2
    assert str(next(SeqIO.parse(out, "fasta")).seq) == "GGGGCCCC"


def test_the_cache_does_not_grow_without_bound(tmp_path, monkeypatch):
    """A batch run walks many genomes; the cache must evict, and evicting must
    close the handle it drops rather than leaking it."""
    genomes = []
    for i in range(6):
        path = tmp_path / f"g{i}.fa"
        path.write_text(f">c1\n{'ACGT' * 100}\n")
        genomes.append(path)
    for path in genomes:
        _extract_window(path, ("c1", 1, 20), tmp_path)
    assert len(search_module._GENOME_INDEX_CACHE) <= search_module._GENOME_INDEX_CACHE_MAX
