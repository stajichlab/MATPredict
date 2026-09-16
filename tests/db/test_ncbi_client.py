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
