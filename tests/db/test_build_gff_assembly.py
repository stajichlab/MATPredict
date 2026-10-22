"""Issue #121: `build_gff_for_record` must keep the CDS feature and protein of an
accession-less gene on an `assembly`-type segment (GCA_/GCF_ accession, fetchable
`seq_region` contig), translating from the recorded exons with the gene's own
codon_start/transl_table -- and must leave accession-bearing assembly records as
they were."""
from __future__ import annotations

from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pytest
import yaml
from Bio import SeqIO
from Bio.Seq import Seq

from MATPredict.db.cli import build_gff_for_record
from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.uniprot_client import UniprotClient
from MATPredict.db.validate import _independent_translation, validate_record

_REPO_DB = Path(__file__).resolve().parents[2] / "db"

# --- synthetic assembly-segment record ---------------------------------------------

_CONTIG_ID = "CONTIG1.1"
_ASSEMBLY = "GCA_000000001.1"
_CONTIG_LEN = 400


def _revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _build_contig() -> str:
    contig = ["C"] * _CONTIG_LEN

    def place(start: int, seq: str) -> None:  # 1-based start
        contig[start - 1 : start - 1 + len(seq)] = list(seq)

    # Plus-strand gene, transl_table 12 (CTG = Ser, not Leu): exons 111-119, 131-139.
    place(111, "ATGCTGAAA")
    place(131, "TTTGGCTAA")
    # Minus-strand gene, codon_start 2: transcript exon 1 (291-299) then exon 2 (260-269).
    # Assembled: G|ATG TGG CA + T CGC AAA TAG -> drop 1 base -> MWHRK*.
    place(291, _revcomp("GATGTGGCA"))
    place(260, _revcomp("TCGCAAATAG"))
    # Accession-bearing single-span gene: 201-212 ATG AAA CCC TAA -> MKP.
    place(201, "ATGAAACCCTAA")
    return "".join(contig)


_CONTIG = _build_contig()
_PLUS_PROTEIN = "MSKFG"
_MINUS_PROTEIN = "MWHRK"
_ACC_PROTEIN = "MKP"


def _transport(requested: list[str]):
    """Fake NCBI E-utilities: nuccore subranges of the contig (strand=2 is the reverse
    complement, as efetch returns it), one protein. Records every requested id."""

    def transport(url: str) -> str:
        params = parse_qs(urlparse(url).query)
        accession = params["id"][0]
        requested.append(accession)
        if params["db"] == ["protein"]:
            assert accession == "ACCPROT1.1"
            return f">{accession}\n{_ACC_PROTEIN}\n"
        assert params["db"] == ["nuccore"]
        if accession != _CONTIG_ID:
            raise AssertionError(f"not a nucleotide record: {accession}")
        start, stop = int(params["seq_start"][0]), int(params["seq_stop"][0])
        seq = _CONTIG[start - 1 : stop]
        if params.get("strand") == ["2"]:
            seq = _revcomp(seq)
        return f">{accession}:{start}-{stop}\n{seq}\n"

    return transport


def _assembly_record(record_id: str) -> dict:
    return {
        "record_id": record_id,
        "mating_type": {"locus_name": "PR"},
        "locus": {"coordinate_provenance": "curator_derived", "core": {"segments": [{
            "segment_index": 0,
            "sequence_source": {"type": "assembly", "accession": _ASSEMBLY, "seq_region": _CONTIG_ID},
            "start": 101, "end": 350,
        }]}},
        "genes": [
            {
                "gene_index": 0, "name": "plusgene", "role": "core_MAT", "present": True,
                "protein_accession": None, "segment_index": 0,
                "start": 111, "end": 139, "strand": "+", "transl_table": 12,
                "exons": [{"start": 111, "end": 119}, {"start": 131, "end": 139}],
            },
            {
                "gene_index": 1, "name": "minusgene", "role": "core_MAT", "present": True,
                "protein_accession": None, "segment_index": 0,
                "start": 260, "end": 299, "strand": "-", "codon_start": 2,
                "exons": [{"start": 291, "end": 299}, {"start": 260, "end": 269}],
            },
            {
                "gene_index": 2, "name": "accgene", "role": "core_MAT", "present": True,
                "protein_accession": "ncbi_protein:ACCPROT1.1", "segment_index": 0,
                "start": 201, "end": 212, "strand": "+",
            },
        ],
    }


def _write(db_root: Path, record: dict) -> Path:
    record_dir = db_root / "Basidiomycota" / "Agaricales" / record["record_id"]
    record_dir.mkdir(parents=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(record, sort_keys=False))
    return record_dir


def _clients(tmp_path: Path, requested: list[str]):
    fetcher = CachedFetcher(cache_dir=tmp_path / "cache", transport=_transport(requested))
    return NcbiClient(email=None, api_key=None, fetcher=fetcher), UniprotClient(fetcher=fetcher)


def _cds_by_gene(gbk_path: Path) -> dict[str, tuple]:
    out = {}
    for seq_record in SeqIO.parse(gbk_path, "genbank"):
        for feature in seq_record.features:
            if feature.type == "CDS":
                out[feature.qualifiers["gene"][0]] = (feature, seq_record)
    return out


def test_accession_less_assembly_genes_keep_cds_and_protein(tmp_path):
    db_root = tmp_path / "db"
    record_dir = _write(db_root, _assembly_record("1_x_PR_B1"))
    requested: list[str] = []
    ncbi, uniprot = _clients(tmp_path, requested)

    build_gff_for_record(db_root, "Basidiomycota", "Agaricales", "1_x_PR_B1", ncbi, uniprot)

    proteins = {r.id.split("|")[2].removeprefix("name="): str(r.seq)
                for r in SeqIO.parse(record_dir / "proteins.faa", "fasta")}
    assert proteins == {"plusgene": _PLUS_PROTEIN, "minusgene": _MINUS_PROTEIN, "accgene": _ACC_PROTEIN}

    cds = _cds_by_gene(record_dir / "locus.gbk")
    assert set(cds) == {"plusgene", "minusgene", "accgene"}
    plus_feature, seq_record = cds["plusgene"]
    assert plus_feature.qualifiers["translation"] == [_PLUS_PROTEIN]
    assert plus_feature.qualifiers["transl_table"] == ["12"]
    assert [(int(p.start), int(p.end), p.strand) for p in plus_feature.location.parts] == [(10, 19, 1), (30, 39, 1)]
    minus_feature, _ = cds["minusgene"]
    assert minus_feature.qualifiers["translation"] == [_MINUS_PROTEIN]
    assert minus_feature.qualifiers["codon_start"] == ["2"]
    assert [(int(p.start), int(p.end), p.strand) for p in minus_feature.location.parts] == [(190, 199, -1), (159, 169, -1)]

    # The written CDS, read back from the written sequence, translates to its own /translation.
    for name, (feature, rec) in cds.items():
        codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
        table = int(feature.qualifiers.get("transl_table", ["1"])[0])
        nt = feature.extract(rec.seq)[codon_start - 1:]
        nt = nt[: len(nt) - len(nt) % 3]
        assert str(nt.translate(table=table, to_stop=True)) == feature.qualifiers["translation"][0], name

    # Every nucleotide fetch went to the seq_region contig, never to the GCA_ accession.
    assert _ASSEMBLY not in requested


def test_independent_translation_assembly_needs_opt_in(tmp_path):
    """`validate_record`'s cross-check keeps the old behaviour (not applicable on an
    assembly segment); only `build_gff_for_record` opts in to the seq_region fetch."""
    requested: list[str] = []
    ncbi, _ = _clients(tmp_path, requested)
    record = _assembly_record("1_x_PR_B1")
    plus, minus = record["genes"][0], record["genes"][1]

    assert _independent_translation(record, plus, ncbi) is None
    assert requested == []
    assert _independent_translation(record, plus, ncbi, resolve_assembly=True) == _PLUS_PROTEIN
    assert _independent_translation(record, minus, ncbi, resolve_assembly=True) == _MINUS_PROTEIN


def test_validate_record_on_assembly_segment_is_unchanged(tmp_path):
    requested: list[str] = []
    ncbi, uniprot = _clients(tmp_path, requested)
    record = {**_assembly_record("1_x_PR_B1"), "taxonomy": {"taxid": 5346}}

    result = validate_record(
        record, ncbi, uniprot,
        taxonomy_runner=lambda cmd, **kw: type("R", (), {"returncode": 0, "stdout": "5346\tk__Fungi\n"})(),
    )

    assert result["sequence_match"]["status"] == "not_applicable"
    assert [g["status"] for g in result["sequence_match"]["per_gene"]] == ["not_applicable"]


# --- regression: committed assembly records regenerate unchanged ---------------------

def _assembly_record_ids() -> list[str]:
    ids = []
    for meta in sorted(_REPO_DB.glob("*/*/*/metadata.yaml")):
        if meta.relative_to(_REPO_DB).parts[0] == "candidates":
            continue
        record = yaml.safe_load(meta.read_text())
        segments = record["locus"].get("core", {}).get("segments", [])
        if any(s["sequence_source"].get("type") == "assembly" for s in segments):
            ids.append("/".join(meta.relative_to(_REPO_DB).parts[:3]))
    return ids


class _OfflineNcbi:
    """Serves each segment's nucleotide sequence from the committed locus.gbk and each
    protein accession's sequence from the committed proteins.faa -- so the committed
    files can be regenerated offline and compared byte for byte."""

    def __init__(self, record: dict, record_dir: Path):
        self._contigs: dict[str, tuple[int, str]] = {}
        gbk_records = list(SeqIO.parse(record_dir / "locus.gbk", "genbank"))
        for segment, seq_record in zip(record["locus"]["core"]["segments"], gbk_records):
            self._contigs[segment["sequence_source"]["seq_region"]] = (segment["start"], str(seq_record.seq))
        by_index = {int(r.id.split("|")[1].removeprefix("gene_index=")): str(r.seq)
                    for r in SeqIO.parse(record_dir / "proteins.faa", "fasta")}
        self._proteins = {
            g["protein_accession"].split(":", 1)[1]: by_index[g["gene_index"]]
            for g in record["genes"] if g.get("protein_accession") and g["gene_index"] in by_index
        }

    def fetch_nucleotide_sequence(self, accession, start, end, strand=None):
        offset, seq = self._contigs[accession]
        sub = seq[start - offset : end - offset + 1]
        return _revcomp(sub) if strand == "-" else sub

    def fetch_protein_sequence(self, accession):
        return self._proteins[accession]


@pytest.mark.parametrize("record_path", [
    r for r in _assembly_record_ids() if not r.endswith("5346_a43-b43-okayama-7_PR_B43")
])
def test_committed_accession_bearing_assembly_records_regenerate_unchanged(tmp_path, record_path):
    phylum, order_or_family, record_id = record_path.split("/")
    src_dir = _REPO_DB / record_path
    record = yaml.safe_load((src_dir / "metadata.yaml").read_text())
    assert all(g.get("protein_accession") for g in record["genes"] if g.get("present", True))

    db_root = tmp_path / "db"
    out_dir = db_root / record_path
    out_dir.mkdir(parents=True)
    (out_dir / "metadata.yaml").write_text((src_dir / "metadata.yaml").read_text())
    (db_root / phylum / "order.yml").write_text((_REPO_DB / phylum / "order.yml").read_text())
    ncbi = _OfflineNcbi(record, src_dir)

    build_gff_for_record(db_root, phylum, order_or_family, record_id, ncbi, uniprot=None)

    for name in ("locus.gbk", "locus.gff3", "proteins.faa"):
        assert (out_dir / name).read_text() == (src_dir / name).read_text(), name


def test_b43_accession_less_genes_regenerate_with_cds_and_protein(tmp_path):
    """The one record exposed by #121: its five accession-less genes (phb3.1, phb1.1,
    phb2.3, phb2.1, phb3.3) keep their CDS and protein, identical to the committed ones."""
    record_path = "Basidiomycota/Agaricales/5346_a43-b43-okayama-7_PR_B43"
    src_dir = _REPO_DB / record_path
    record = yaml.safe_load((src_dir / "metadata.yaml").read_text())
    db_root = tmp_path / "db"
    out_dir = db_root / record_path
    out_dir.mkdir(parents=True)
    (out_dir / "metadata.yaml").write_text((src_dir / "metadata.yaml").read_text())
    ncbi = _OfflineNcbi(record, src_dir)

    build_gff_for_record(db_root, *record_path.split("/"), ncbi, uniprot=None)

    assert (out_dir / "proteins.faa").read_text() == (src_dir / "proteins.faa").read_text()
    cds = _cds_by_gene(out_dir / "locus.gbk")
    assert len(cds) == 12
    for name in ("phb3.1", "phb1.1", "phb2.3", "phb2.1", "phb3.3"):
        assert name in cds
