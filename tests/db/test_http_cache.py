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


def test_the_cache_key_ignores_the_caller_identity(tmp_path):
    """E-mail, API key and tool name identify the caller, not the request, so
    one cache serves every user's configuration."""
    calls = []
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=lambda u: calls.append(u) or "body")
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=4837&retmode=xml"
    fetcher.get(base + "&tool=MATPredict&email=a@example.org")
    fetcher.get(base + "&tool=MATPredict&email=b@example.org&api_key=XYZ")
    fetcher.get(base + "&tool=MATPredict")
    assert len(calls) == 1


def test_cache_key_keeps_every_other_byte():
    from MATPredict.db.http_cache import cache_key
    assert cache_key("https://x.org/a?db=n&id=A,B:1&email=e@x&strand=2") == "https://x.org/a?db=n&id=A,B:1&strand=2"
    assert cache_key("https://x.org/a") == "https://x.org/a"
    assert cache_key("https://x.org/a?email=e") == "https://x.org/a"


def test_an_old_cache_entry_is_found_and_copied_forward(tmp_path):
    """Caches written before 2026-10-03 are keyed on the full URL, which ended
    with the fixed default e-mail. They must still answer, with no request."""
    import hashlib
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=4837&retmode=xml"
    old = tmp_path / (hashlib.sha256((base + "&email=jason.stajich@ucr.edu").encode()).hexdigest() + ".txt")
    old.write_text("cached body")
    calls = []
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=lambda u: calls.append(u) or "fresh")
    assert fetcher.get(base + "&tool=MATPredict&email=someone@example.org") == "cached body"
    assert calls == []
    assert fetcher._cache_path(base).read_text() == "cached body"


def test_the_old_key_is_not_tried_for_other_hosts(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=lambda u: "x")
    assert fetcher._legacy_path("https://rest.uniprot.org/uniprotkb/P12345.fasta") is None
