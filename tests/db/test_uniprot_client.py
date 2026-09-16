from __future__ import annotations

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.uniprot_client import UniprotClient

ENTRY_JSON = """{"primaryAccession": "P12345", "entryType": "reviewed"}"""
FASTA = ">sp|P12345|SEXP_PHYBL Sex pheromone\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n"


def _fake_transport(responses):
    def transport(url: str) -> str:
        for key, body in responses.items():
            if key in url:
                return body
        raise AssertionError(f"unexpected url: {url}")
    return transport


def test_resolve_accession_found(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"uniprotkb/P12345.json": ENTRY_JSON}))
    client = UniprotClient(fetcher=fetcher)
    status = client.resolve_accession("P12345")
    assert status.resolved is True
    assert status.resolved_version == "P12345"


def test_fetch_protein_sequence(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"uniprotkb/P12345.fasta": FASTA}))
    client = UniprotClient(fetcher=fetcher)
    seq = client.fetch_protein_sequence("P12345")
    assert seq.startswith("MKTAYIAKQRQ")
