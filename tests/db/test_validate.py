from __future__ import annotations

from types import SimpleNamespace

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.uniprot_client import UniprotClient
from MATPredict.db.validate import validate_record

RECORD = {
    "taxonomy": {"taxid": 4837},
    "locus": {
        "coordinate_provenance": "published_explicit",
        "core": {"segments": [{"segment_index": 0, "sequence_source": {"type": "assembly", "accession": "GCA_000315115.1"}}]},
    },
    "genes": [
        {"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1", "present": True},
    ],
}

ESUMMARY_LIVE = '{"result": {"uids": ["1"], "1": {"accessionversion": "GCA_000315115.1", "status": "live"}}}'
EFETCH_FASTA = ">AAB12345.1\nMKTAYIAKQRQISFVKSHFSRQ\n"


def _fake_transport(url: str) -> str:
    if "esummary" in url:
        return ESUMMARY_LIVE
    if "efetch" in url:
        return EFETCH_FASTA
    raise AssertionError(url)


def _fake_taxonomy_runner(cmd, **kwargs):
    return SimpleNamespace(returncode=0, stdout="4837\tk__Fungi;p__Mucoromycota;...;s__Phycomyces_blakesleeanus\n")


def test_validate_record_all_pass(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport)
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    uniprot = UniprotClient(fetcher=fetcher)

    result = validate_record(
        record={**RECORD, "genes": [{**RECORD["genes"][0], "protein_accession": "ncbi_protein:AAB12345.1"}]},
        ncbi=ncbi,
        uniprot=uniprot,
        taxonomy_runner=_fake_taxonomy_runner,
    )

    assert result["accession_resolved"] is True
    assert result["accession_resolved_version"] == "GCA_000315115.1"
    assert result["sequence_match"]["status"] in {"pass", "warn"}
    assert result["taxonomy_current"] is True


def test_validate_record_skips_coordinate_checks_when_not_available():
    record = {
        "taxonomy": {"taxid": 4837},
        "locus": {"coordinate_provenance": "not_available"},
        "genes": [],
    }
    result = validate_record(record=record, ncbi=None, uniprot=None, taxonomy_runner=_fake_taxonomy_runner)
    assert result["accession_resolved"] is None
    assert result["taxonomy_current"] is True
