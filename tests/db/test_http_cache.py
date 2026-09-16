from __future__ import annotations

from MATPredict.db.http_cache import CachedFetcher


def test_cache_avoids_second_transport_call(tmp_path):
    calls = []

    def transport(url: str) -> str:
        calls.append(url)
        return f"response for {url}"

    fetcher = CachedFetcher(cache_dir=tmp_path, transport=transport)
    first = fetcher.get("https://example.org/a")
    second = fetcher.get("https://example.org/a")

    assert first == second == "response for https://example.org/a"
    assert calls == ["https://example.org/a"]  # only called once


def test_cache_distinguishes_urls(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=lambda url: url)
    assert fetcher.get("https://example.org/a") != fetcher.get("https://example.org/b")
