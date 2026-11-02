"""`--phylum auto`: the decision table, the GmG record reader, the index cache, and the routing wired into `_cmd_detect`."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

import MATPredict.detect.cli as detect_cli
from MATPredict.detect import gmg_guess
from MATPredict.detect.family_registry import Family, FamilyKey, _gene_alias_map, route
from MATPredict.detect.gmg_guess import GmgGuess, decide, guess_from_record, newer_index, resolve_index

AVAILABLE = ["Ascomycota", "Basidiomycota", "Mucoromycota"]


def _call(taxon, taxid=None):
    return {"taxon": taxon, "taxid": taxid}


def _rec(phylum="Ascomycota", flags=(), family=("Nectriaceae", 110618), order=("Hypocreales", 5125), hint=None, status="ok"):
    calls = {"phylum": _call(phylum, 4890 if phylum else None), "class": _call(None), "order": _call(*order),
             "family": _call(*family), "genus": _call("Fusarium", 5506), "species": _call(None)}
    rec = {"status": status, "qc_flags": list(flags), "guess": {"calls": calls}}
    if hint is not None:
        rec["routing_hint"] = hint
    return rec


# ---- the record reader

def test_guess_uses_the_routing_hint_when_the_record_has_one():
    hint = {"usable": True, "reason": None, "phylum": "Ascomycota", "phylum_taxid": 4890, "deepest_rank": "genus"}
    g = guess_from_record(_rec(hint=hint))
    assert (g.usable, g.phylum, g.deepest_rank) == (True, "Ascomycota", "genus")


def test_routing_taxid_is_the_family_then_the_order_never_the_genus():
    assert guess_from_record(_rec()).routing_taxid == 110618
    g = guess_from_record(_rec(family=(None, None)))
    assert (g.routing_taxid, g.routing_taxid_rank) == (5125, "order")
    g = guess_from_record(_rec(family=(None, None), order=(None, None)))
    assert g.routing_taxid is None


@pytest.mark.parametrize("flags,status,phylum,reason", [
    ((), "ok", "Ascomycota", None),
    (("weak_phylum",), "low_confidence", "Ascomycota", "weak_phylum"),
    (("not_fungal",), "no_match", "Chordata", "not_fungal"),
    ((), "no_match", None, "no_phylum_call"),
    ((), "error", None, "error"),
])
def test_a_record_without_a_hint_gets_the_same_rule(flags, status, phylum, reason):
    g = guess_from_record(_rec(phylum=phylum, flags=flags, status=status))
    assert g.reason == reason and g.usable == (reason is None)


# ---- the decision

def test_a_usable_curated_phylum_is_routed_with_its_taxid():
    d = decide(GmgGuess(usable=True, phylum="Mucoromycota", routing_taxid=4827), AVAILABLE)
    assert (d.phylum, d.taxid, d.weak) == ("Mucoromycota", 4827, False)


@pytest.mark.parametrize("guess,needle", [
    (GmgGuess(reason="error"), "could not read"),
    (GmgGuess(reason="no_phylum_call"), "no phylum call"),
    (GmgGuess(reason="not_fungal", phylum="Chordata"), "not a fungus"),
    (GmgGuess(reason="weak_phylum", phylum="Ascomycota"), "auto-weak"),
    (GmgGuess(usable=True, phylum="Chytridiomycota"), "no curated MAT locus family"),
    (GmgGuess(reason="something_new", phylum="Ascomycota"), "not usable"),
])
def test_anything_doubtful_is_not_searched_and_says_why(guess, needle):
    d = decide(guess, AVAILABLE)
    assert d.phylum is None and needle in d.note


def test_a_weak_call_is_searched_only_with_auto_weak_and_is_marked():
    g = GmgGuess(reason="weak_phylum", phylum="Ascomycota", routing_taxid=5125)
    assert decide(g, AVAILABLE).phylum is None
    d = decide(g, AVAILABLE, accept_weak=True)
    assert d.phylum == "Ascomycota" and d.weak and "weak" in d.note
    # auto-weak never accepts a non-fungal call or a phylum without families
    assert decide(GmgGuess(reason="not_fungal", phylum="Chordata"), AVAILABLE, accept_weak=True).phylum is None
    assert decide(GmgGuess(reason="weak_phylum", phylum="Chytridiomycota"), AVAILABLE, accept_weak=True).phylum is None


# ---- the index

def test_newer_index_only_for_a_newer_representatives_index():
    manifest = [{"file": "bfd-reps-2026.10.5-s20000.gmgi"}, {"file": "bfd-reps-2026.11-s20000.gmgi"},
                {"file": "bfd-2027.1-s10000.gmgi"}]
    assert newer_index("bfd-reps-2026.10.5-s20000.gmgi", manifest) == "bfd-reps-2026.11-s20000.gmgi"
    assert newer_index("bfd-reps-2026.11-s20000.gmgi", manifest) is None
    assert newer_index("bfd-reps-2026.10.5-s20000.gmgi", None) is None


def test_resolve_index_returns_an_existing_path_untouched(tmp_path):
    f = tmp_path / "x.gmgi"
    f.write_bytes(b"x")
    assert resolve_index(str(f)) == (f.resolve(), str(f))


def test_resolve_index_reports_a_missing_path(tmp_path):
    with pytest.raises(gmg_guess.GmgUnavailable, match="not found"):
        resolve_index(str(tmp_path / "nope" / "x.gmgi"))


def test_resolve_index_downloads_once_into_the_cache_and_checks_the_checksum(tmp_path, monkeypatch):
    calls = []

    def fake_download(url, dest):
        calls.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"index bytes")

    monkeypatch.setattr(gmg_guess, "_download", fake_download)
    monkeypatch.setattr(gmg_guess, "_manifest", lambda cache: None)
    monkeypatch.setattr(gmg_guess, "_release_sha256", lambda name: gmg_guess._sha256(tmp_path / name))
    name = "bfd-reps-2026.10.5-s20000.gmgi"
    # a published checksum that does not match the download is an error and the file is removed
    monkeypatch.setattr(gmg_guess, "_release_sha256", lambda n: "0" * 64)
    with pytest.raises(gmg_guess.GmgUnavailable, match="sha256"):
        resolve_index(name, cache=tmp_path)
    assert not (tmp_path / name).exists()
    import hashlib
    monkeypatch.setattr(gmg_guess, "_release_sha256", lambda n: hashlib.sha256(b"index bytes").hexdigest())
    path, spec = resolve_index(name, cache=tmp_path)
    assert path == tmp_path / name and spec == name
    resolve_index(name, cache=tmp_path)           # cached: no second download
    assert len(calls) == 2 and calls[0] == gmg_guess.GMG_INDEX_SITE + name


def test_resolve_index_default_is_the_pinned_index_or_the_environment(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(gmg_guess, "_manifest", lambda cache: None)
    monkeypatch.setattr(gmg_guess, "_download", lambda url, dest: (seen.append(dest.name), dest.write_bytes(b"x")))
    monkeypatch.setattr(gmg_guess, "_release_sha256", lambda n: None)
    monkeypatch.delenv("MATPREDICT_GMG_INDEX", raising=False)
    resolve_index(None, cache=tmp_path)
    monkeypatch.setenv("MATPREDICT_GMG_INDEX", "bfd-2026.10.5-s10000.gmgi")
    resolve_index(None, cache=tmp_path)
    assert seen == [gmg_guess.DEFAULT_GMG_INDEX, "bfd-2026.10.5-s10000.gmgi"]


def test_cache_dir_precedence(monkeypatch, tmp_path):
    for v in ("MATPREDICT_GMG_CACHE", "MATPREDICT_CACHE_DIR", "XDG_CACHE_HOME"):
        monkeypatch.delenv(v, raising=False)
    assert gmg_guess.gmg_cache_dir() == Path("~/.cache").expanduser() / "matpredict" / "gmg"
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "x"))
    assert gmg_guess.gmg_cache_dir() == tmp_path / "x" / "matpredict" / "gmg"
    monkeypatch.setenv("MATPREDICT_CACHE_DIR", str(tmp_path / "c"))
    assert gmg_guess.gmg_cache_dir() == tmp_path / "c" / "gmg"
    monkeypatch.setenv("MATPREDICT_GMG_CACHE", str(tmp_path / "g"))
    assert gmg_guess.gmg_cache_dir() == tmp_path / "g"


# ---- routing wired into the command

def _fam(phylum, name, scope, sub=None):
    genes = [{"name": "G", "role": "core_MAT"}]
    return Family(key=FamilyKey(phylum, name), vocabulary_type="enum", idiomorph_values=["a"], idiomorph_pattern=None,
                  genes=genes, taxonomic_scope=scope, gene_aliases=_gene_alias_map(genes), fallback_subphylum=sub)


@pytest.fixture
def wired(monkeypatch):
    """`_route_with_gmg` over two small fake phyla; GmG and the taxonomy lookup are faked."""
    fams = [_fam("Ascomycota", "MAT", [716546]), _fam("Basidiomycota", "HD", [5204])]
    monkeypatch.setattr(detect_cli, "load_all_families", lambda root: fams)
    monkeypatch.setattr(detect_cli, "available_phyla", lambda root: ["Ascomycota", "Basidiomycota"])
    lineages = {5125: [716546, 4890], 999: [5204]}
    real = detect_cli.route
    monkeypatch.setattr(detect_cli, "route", lambda taxid, families, **kw: real(
        taxid, families, lineage_taxids_resolver=lambda t: lineages[t], phylum_name_resolver=lambda t: None, **kw))
    state = SimpleNamespace(guess=None, ran=0)

    def fake_guess(genome, index):
        state.ran += 1
        return state.guess
    monkeypatch.setattr(detect_cli, "guess_phylum", fake_guess)
    cfg = SimpleNamespace(db_root=Path("."))

    def go(phylum="auto", taxid=None, exhaustive=False):
        args = SimpleNamespace(genome="g.fa", taxid=taxid, phylum=phylum, gmg_index=None, exhaustive=exhaustive,
                               scoped_roster_policy=None)
        monkeypatch.setattr(detect_cli, "_scoped_roster_policy_from_args", lambda a: "off")
        return args, detect_cli._route_with_gmg(args, cfg)
    state.go = go
    return state


def test_a_usable_call_with_a_taxid_routes_by_lineage(wired):
    wired.guess = GmgGuess(usable=True, phylum="Ascomycota", routing_taxid=5125, routing_taxid_rank="order")
    args, routing = wired.go()
    assert routing.routing_mode == "lineage" and [f.key.phylum for f in routing.families] == ["Ascomycota"]
    assert args.phylum == "Ascomycota" and args._phylum_guess["action"] == "routed"
    assert "order taxid 5125" in args._phylum_guess["routed_by"]


def test_the_genetic_code_is_read_from_the_gmg_taxid(wired, monkeypatch):
    wired.guess = GmgGuess(usable=True, phylum="Ascomycota", routing_taxid=5125, routing_taxid_rank="order")
    monkeypatch.setattr(detect_cli, "default_genetic_code", lambda t: 12 if t == 5125 else 1)
    args, routing = wired.go()
    assert detect_cli._resolve_genetic_code(args) == 12
    assert args._phylum_guess["genetic_code"] == 12 and "order taxid 5125" in args._phylum_guess["genetic_code_from"]


def test_an_explicit_genetic_code_wins_and_a_failed_lookup_is_recorded(wired, monkeypatch):
    wired.guess = GmgGuess(usable=True, phylum="Ascomycota", routing_taxid=5125, routing_taxid_rank="order")

    def boom(t):
        raise OSError("no network")
    monkeypatch.setattr(detect_cli, "default_genetic_code", boom)
    args, _ = wired.go()
    assert detect_cli._resolve_genetic_code(args) is None and "used table 1" in args._phylum_guess["genetic_code_error"]
    monkeypatch.setattr(detect_cli, "default_genetic_code", lambda t: 12)
    args = SimpleNamespace(genome="g.fa", taxid=None, phylum="auto", gmg_index=None, exhaustive=False, genetic_code=4)
    detect_cli._route_with_gmg(args, SimpleNamespace(db_root=Path(".")))
    assert detect_cli._resolve_genetic_code(args) == 4 and "genetic_code" not in args._phylum_guess


def test_a_taxid_that_lands_outside_the_called_phylum_is_ignored(wired):
    wired.guess = GmgGuess(usable=True, phylum="Ascomycota", routing_taxid=999, routing_taxid_rank="family")
    args, routing = wired.go()
    assert routing.routing_mode == "explicit_phylum" and [f.key.phylum for f in routing.families] == ["Ascomycota"]


def test_a_usable_call_without_a_taxid_routes_by_phylum(wired):
    wired.guess = GmgGuess(usable=True, phylum="Basidiomycota")
    args, routing = wired.go()
    assert routing.routing_mode == "explicit_phylum" and args.phylum == "Basidiomycota"


def test_a_weak_call_is_not_searched_unless_auto_weak(wired):
    wired.guess = GmgGuess(reason="weak_phylum", phylum="Ascomycota")
    args, routing = wired.go()
    assert routing.routing_mode == "not_searched" and routing.families == [] and args.phylum is None
    assert "weak" in detect_cli._not_searched_reason(args, routing)
    args, routing = wired.go(phylum="auto-weak")
    assert routing.routing_mode == "explicit_phylum" and args._phylum_guess["weak_call_accepted"] is True


def test_not_searched_by_gmg_can_still_be_forced_exhaustive(wired):
    wired.guess = GmgGuess(reason="not_fungal", phylum="Chordata")
    args, routing = wired.go(exhaustive=True)
    assert routing.routing_mode == "exhaustive" and len(routing.families) == 2


def test_an_explicit_taxid_wins_and_gmg_is_not_run(wired):
    args, routing = wired.go(taxid=5125)
    assert wired.ran == 0 and routing.routing_mode == "lineage"
    assert args._phylum_guess["action"] == "skipped" and args.phylum is None


def test_parser_offers_auto_and_the_gmg_index(monkeypatch, tmp_path):
    from MATPredict.__main__ import build_parser
    (tmp_path / "Ascomycota").mkdir()
    (tmp_path / "Ascomycota" / "order.yml").write_text("phylum: Ascomycota\nloci: []\n")
    monkeypatch.setenv("MATPREDICT_DB_ROOT", str(tmp_path))
    a = build_parser().parse_args(["detect", "--genome", "g.fa", "--phylum", "auto-weak", "--gmg-index", "latest"])
    assert (a.phylum, a.gmg_index) == ("auto-weak", "latest")


def test_an_unrelated_gmg_package_on_pypi_is_not_mistaken_for_ours(monkeypatch):
    import sys
    import types
    monkeypatch.setitem(sys.modules, "gmg", types.ModuleType("gmg"))      # no `Index`: not GmG's bindings
    monkeypatch.setattr(gmg_guess.shutil, "which", lambda name: None)
    with pytest.raises(gmg_guess.GmgUnavailable, match="not `pip install gmg`"):
        gmg_guess.classify_genome("g.fa", Path("x.gmgi"))
