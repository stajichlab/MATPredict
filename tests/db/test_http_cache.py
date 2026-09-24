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


def test_an_empty_cache_file_is_a_miss_not_an_answer(tmp_path):
    """A reader that raced a writer could find a truncated file. Returning ''
    as the NCBI answer made the taxonomy parse fail, and routing silently
    fell back to exhaustive on 128 concurrent pilot runs."""
    calls = []
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=lambda u: calls.append(u) or "body")
    fetcher._cache_path("https://example.org/a").parent.mkdir(parents=True, exist_ok=True)
    fetcher._cache_path("https://example.org/a").write_text("")
    assert fetcher.get("https://example.org/a") == "body"
    assert calls == ["https://example.org/a"]


def test_a_write_leaves_no_temporary_file(tmp_path):
    """The body is written to a temporary name and renamed into place, so a
    concurrent reader sees either nothing or the whole response."""
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=lambda u: "body")
    fetcher.get("https://example.org/a")
    assert [p.name for p in tmp_path.iterdir()] == [fetcher._cache_path("https://example.org/a").name]
