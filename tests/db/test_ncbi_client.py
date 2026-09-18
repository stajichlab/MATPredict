from __future__ import annotations

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient

ESUMMARY_LIVE = """{"result": {"uids": ["1"], "1": {"accessionversion": "GCA_000315115.1", "status": "live"}}}"""
ESUMMARY_SUPPRESSED = """{"result": {"uids": ["1"], "1": {"accessionversion": "GCA_999999999.1", "status": "suppressed"}}}"""
EFETCH_FASTA = ">AAB12345.1 sexP protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n"


def _fake_transport(responses):
    def transport(url: str) -> str:
        for key, body in responses.items():
            if key in url:
                return body
        raise AssertionError(f"unexpected url: {url}")
    return transport


def test_resolve_accession_live(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"esummary": ESUMMARY_LIVE}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    status = client.resolve_accession("GCA_000315115.1")
    assert status.resolved is True
    assert status.suppressed is False
    assert status.resolved_version == "GCA_000315115.1"


def test_resolve_accession_suppressed(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"esummary": ESUMMARY_SUPPRESSED}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    status = client.resolve_accession("GCA_999999999.1")
    assert status.resolved is False
    assert status.suppressed is True


def test_fetch_protein_sequence(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"efetch": EFETCH_FASTA}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    seq = client.fetch_protein_sequence("AAB12345.1")
    assert seq.startswith("MKTAYIAKQRQ")
    assert "\n" not in seq


def test_fetch_nucleotide_sequence_plus_strand(tmp_path):
    fetcher = CachedFetcher(
        cache_dir=tmp_path,
        transport=_fake_transport({"db=nuccore": ">EU009461.1:100-109\nACGTACGTAC\n"}),
    )
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    seq = client.fetch_nucleotide_sequence("EU009461.1", 100, 109, "+")
    assert seq == "ACGTACGTAC"


EFETCH_TAXONOMY_XML = """<?xml version="1.0"?>
<TaxaSet>
<Taxon>
<TaxId>5270</TaxId>
<ScientificName>Ustilago maydis</ScientificName>
<LineageEx>
<Taxon><TaxId>4751</TaxId><ScientificName>Fungi</ScientificName><Rank>kingdom</Rank></Taxon>
<Taxon><TaxId>5204</TaxId><ScientificName>Basidiomycota</ScientificName><Rank>phylum</Rank></Taxon>
<Taxon><TaxId>5157</TaxId><ScientificName>Ustilaginomycotina</ScientificName><Rank>subphylum</Rank></Taxon>
<Taxon><TaxId>5259</TaxId><ScientificName>Ustilaginomycetes</ScientificName><Rank>class</Rank></Taxon>
</LineageEx>
</Taxon>
</TaxaSet>
"""


def test_fetch_taxonomy_lineage_parses_lineage_ex_taxids(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"db=taxonomy": EFETCH_TAXONOMY_XML}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    lineage = client.fetch_taxonomy_lineage(5270)
    assert lineage == [4751, 5204, 5157, 5259]


def test_fetch_nucleotide_sequence_minus_strand_requests_strand2(tmp_path):
    captured_urls = []

    def transport(url: str) -> str:
        captured_urls.append(url)
        return ">EU009461.1:100-109 c\nGTACGTACGT\n"

    fetcher = CachedFetcher(cache_dir=tmp_path, transport=transport)
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    seq = client.fetch_nucleotide_sequence("EU009461.1", 100, 109, "-")
    assert seq == "GTACGTACGT"
    assert "strand=2" in captured_urls[0]
