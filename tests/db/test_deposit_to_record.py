"""scripts/deposit_to_record.py: deposit -> core-only record, offline against a cached real deposit.

The fixture cache holds the six NCBI responses for Ascochyta rabiei ABC70320.1
(GenPept, protein FASTA, nuccore summary and GenBank window of DQ341313.1,
taxonomy 5454, PubMed summary), exactly as the generator cached them on
2026-10-22, so the test needs no network.
"""
import importlib.util
import shutil
import sys
from pathlib import Path

import pytest
import yaml

from MATPredict.db.identifiers import slugify_strain
from MATPredict.db.schema import (
    load_metadata_schema,
    validate_gene_status,
    validate_gene_vocabulary,
    validate_idiomorphs,
    validate_metadata,
)

_ROOT = Path(__file__).resolve().parents[2]
_CACHE = Path(__file__).parent / "data" / "deposit_to_record_cache"
_spec = importlib.util.spec_from_file_location("deposit_to_record", _ROOT / "scripts" / "deposit_to_record.py")
gen = importlib.util.module_from_spec(_spec)
sys.modules["deposit_to_record"] = gen  # dataclasses resolve annotations through sys.modules
_spec.loader.exec_module(gen)


def _empty_db(tmp: Path) -> Path:
    """A db root holding only the real Ascomycota order.yml (no records to collide with)."""
    root = tmp / "db"
    (root / "Ascomycota").mkdir(parents=True, exist_ok=True)
    shutil.copy(_ROOT / "db" / "Ascomycota" / "order.yml", root / "Ascomycota" / "order.yml")
    return root


def _run(cache_dir: Path, db_root: Path, **kw):
    fetcher = gen.Fetcher(cache_dir, offline=True)
    kw = {"mating_system": "heterothallic", "tier_policy": "cited2-uncited3", **kw}
    return gen.build("ABC70320.1", fetcher, db_root, gen.load_db_index(db_root), date="2026-10-22", **kw)


def _edit_genbank(cache_dir: Path, old: str, new: str) -> None:
    gbk = _file_starting(cache_dir, "LOCUS       DQ341313")
    text = gbk.read_text()
    assert text.count(old) == 1, old
    gbk.write_text(text.replace(old, new))


def _write_record(db: Path, rec: dict) -> None:
    d = db / "Ascomycota" / "Pleosporales" / rec["record_id"]
    d.mkdir(parents=True)
    (d / "metadata.yaml").write_text(yaml.safe_dump(rec, sort_keys=False))


def _copy_cache(tmp_path: Path) -> Path:
    dst = tmp_path / "cache"
    shutil.copytree(_CACHE, dst)
    return dst


def _file_starting(cache_dir: Path, prefix: str) -> Path:
    return next(p for p in cache_dir.iterdir() if p.read_text().startswith(prefix))


def test_real_deposit_gives_a_valid_confirmed_record(tmp_path):
    res = _run(_CACHE, _empty_db(tmp_path))
    assert res.verdict == gen.CONFIRMED, res.reasons("needs_human") + res.reasons("caveat")
    rec = res.record
    assert rec["record_id"] == "5454_atcc-76502_MAT_MAT1-1"
    gene = rec["genes"][0]
    # DQ341313.1 CDS join(1315..1558,1611..2467): transcript order, stop codon inside the last exon
    assert gene["exons"] == [{"start": 1315, "end": 1558}, {"start": 1611, "end": 2467}]
    assert (gene["start"], gene["end"], gene["strand"]) == (1315, 2467, "+")
    assert gene["protein_accession"] == "ncbi_protein:ABC70320.1"
    assert rec["locus"]["core"]["segments"][0]["sequence_source"]["accession"] == "DQ341313.1"
    assert rec["taxonomy"]["lineage"].startswith("k__Fungi;p__Ascomycota;subphylum__Pezizomycotina;")
    assert rec["taxonomy"]["lineage"].endswith(";s__Ascochyta_rabiei")
    assert rec["organism"]["strain"] == {"name": "ATCC 76502", "known": True,
                                         "culture_collection_ids": [], "differs_from_sequenced": False}
    assert rec["evidence"]["locus_existence"]["citations"][0]["pmid"] == "12781674"
    assert rec["mating_type"] == {"locus_name": "MAT", "idiomorphs": ["MAT1-1"], "system": "heterothallic"}
    assert rec["validation"]["status"] == "needs_review"
    assert "CURATOR-SUPPLIED" in rec["curation"]["notes"]
    # ruling 13: core-only record fields
    assert rec["organism"]["species"] == "Ascochyta rabiei"
    assert rec["locus"]["core"]["completeness"] == "partial"
    assert rec["locus"]["excluded_from_coordinate_benchmark"] is True
    seg = rec["locus"]["core"]["segments"][0]
    assert (seg["start"], seg["end"]) == (1315, 2467)  # segment = CDS span
    # ruling 2: cited deposit -> tier 2
    assert {rec["evidence"][k]["tier"] for k in rec["evidence"]} == {2}
    # ruling 6: the annotation source is stated explicitly
    assert "ANNOTATION SOURCE: locus-specific deposit" in rec["curation"]["notes"]
    # the record passes the same checks `curate-db propose/accept/validate` apply
    order_doc = yaml.safe_load((_ROOT / "db" / "Ascomycota" / "order.yml").read_text())
    assert validate_metadata(rec) == []
    assert validate_idiomorphs(rec, order_doc) == []
    assert validate_gene_vocabulary(rec, order_doc) == []
    assert validate_gene_status(rec) == []


def test_output_is_deterministic(tmp_path):
    db = _empty_db(tmp_path)
    a = yaml.safe_dump(_run(_CACHE, db).record, sort_keys=False)
    b = yaml.safe_dump(_run(_CACHE, db).record, sort_keys=False)
    assert a == b


def test_an_already_curated_deposit_needs_a_human(tmp_path):
    db = _empty_db(tmp_path)
    rec = _run(_CACHE, db).record
    d = db / "Ascomycota" / "Pleosporales" / rec["record_id"]
    d.mkdir(parents=True)
    (d / "metadata.yaml").write_text(yaml.safe_dump(rec, sort_keys=False))
    res = _run(_CACHE, db)
    assert res.verdict == gen.NEEDS_HUMAN
    assert any(r.startswith(("record_id:", "duplicate:")) for r in res.reasons("needs_human"))


def test_a_protein_that_the_exons_do_not_encode_needs_a_human(tmp_path):
    cache = _copy_cache(tmp_path)
    fasta = _file_starting(cache, ">ABC70320.1")
    header, seq = fasta.read_text().split("\n", 1)
    seq = seq.replace("\n", "")
    mutated = seq[:20] + ("A" if seq[20] != "A" else "G") + seq[21:]
    fasta.write_text(header + "\n" + mutated + "\n")
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.NEEDS_HUMAN
    assert any(r.startswith("translation:") for r in res.reasons("needs_human"))


def test_a_non_canonical_intron_needs_a_human(tmp_path):
    cache = _copy_cache(tmp_path)
    gbk = _file_starting(cache, "LOCUS       DQ341313")
    lines = gbk.read_text().split("\n")
    # the intron is 1559..1610; its donor G is position 1559 (ORIGIN line 1501, block 6, base 9)
    i = next(n for n, line in enumerate(lines) if line.strip().startswith("1501 "))
    fields = lines[i].split()
    block = fields[6]
    assert block[8:10] == "gt"
    fields[6] = block[:8] + "a" + block[9:]
    lines[i] = lines[i][:10] + " ".join(fields[1:])
    gbk.write_text("\n".join(lines))
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.NEEDS_HUMAN
    assert any("AT..AG" in r for r in res.reasons("needs_human"))


def test_slug_rule_matches_the_db_identifier_rule():
    for name in ["ATCC 76502", "CBS 101121", "UVas1-58", "CMW 39103", "f. sp. x/y", "MLAX_297"]:
        assert gen.slugify_strain(name) == slugify_strain(name, True, set())


def test_never_sends_an_email(tmp_path, monkeypatch):
    seen = []

    class Boom(Exception):
        pass

    def fake_urlopen(url, timeout=0):
        seen.append(url)
        raise Boom

    monkeypatch.setattr(gen.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(gen.time, "sleep", lambda s: None)
    fetcher = gen.Fetcher(tmp_path)
    with pytest.raises(RuntimeError):
        fetcher.get("efetch.fcgi", db="protein", id="X1.1", rettype="fasta", retmode="text")
    assert seen and all("email" not in u for u in seen) and all("tool=" in u for u in seen)


# --------------------------------------------------------------------------- curator rulings 2026-10-22


def test_mating_system_defaults_to_unknown_and_is_never_guessed(tmp_path):
    fetcher = gen.Fetcher(_CACHE, offline=True)
    db = _empty_db(tmp_path)
    res = gen.build("ABC70320.1", fetcher, db, gen.load_db_index(db), date="2026-10-22")
    assert res.verdict == gen.CONFIRMED
    assert res.record["mating_type"]["system"] == "unknown"
    assert "mating_type.system=unknown" in res.record["curation"]["notes"]
    assert "unknown" in gen.MATING_SYSTEMS
    # Schema PR #127 (branch system-unknown) adds `unknown` to the enum. Read the enum from the
    # schema so this test holds before and after it lands: once accepted, the record is valid;
    # until then the system value is the record's only schema error.
    enum = load_metadata_schema()["properties"]["mating_type"]["properties"]["system"]["enum"]
    errors = validate_metadata(res.record)
    if "unknown" in enum:
        assert errors == []
    else:
        assert len(errors) == 1 and "unknown" in errors[0], errors


def test_cli_defaults_mating_system_unknown_and_tier_policy(tmp_path):
    repo = tmp_path
    _empty_db(repo)
    report = tmp_path / "rep" / "r"
    assert gen.main(["--repo", str(repo), "--protein", "ABC70320.1", "--date", "2026-10-22",
                     "--cache-dir", str(_CACHE), "--offline", "--out", str(tmp_path / "out"),
                     "--report", str(report)]) == 0
    rec = yaml.safe_load(next((tmp_path / "out").glob("**/metadata.yaml")).read_text())
    assert rec["mating_type"]["system"] == "unknown"
    assert "cited2-uncited3 (tier 2)" in rec["curation"]["notes"]
    assert (tmp_path / "rep" / "r.conflicts.tsv").read_text().split("\n")[0] == "\t".join(gen.CONFLICT_COLUMNS)


def test_builtin_hold_list_holds_helotiales_mat1_1_1_only():
    holds = gen.load_holds(None)
    held = gen.matching_holds(holds, "QRL06102.1", "Helotiales", "MAT1-1-1")  # Diplocarpon coronariae
    assert [h["reason"] for h in held] == [
        "held: Helotiales MAT1-1-1 until the idiomorph vote ignores short opposite-idiomorph fragments"]
    assert gen.matching_holds(holds, "QRL06106.1", "Helotiales", "MAT1-2-1") == []  # its MAT1-2-1: not held
    assert gen.matching_holds(holds, "ABC70320.1", "Pleosporales", "MAT1-1-1") == []
    assert gen.matching_holds(gen.load_holds(None, builtin=False), "QRL06102.1", "Helotiales", "MAT1-1-1") == []


def test_hold_rules_make_a_deposit_need_a_human(tmp_path):
    db = _empty_db(tmp_path)
    by_order = [{"protein": "", "order": "Pleosporales", "gene": "MAT1-1-1", "reason": "held: test order rule"}]
    res = _run(_CACHE, db, holds=by_order)
    assert res.verdict == gen.NEEDS_HUMAN and res.reasons("needs_human") == ["hold: held: test order rule"]
    hold_file = tmp_path / "holds.tsv"
    hold_file.write_text("# extra holds\nprotein\torder\tgene\treason\nABC70320\t\t\t\n")
    res = _run(_CACHE, db, holds=gen.load_holds(hold_file))
    assert res.verdict == gen.NEEDS_HUMAN and res.reasons("needs_human") == ["hold: held: listed in holds.tsv"]
    hold_file.write_text("protein\torder\tgene\treason\n\t\t\tholds everything\n")
    with pytest.raises(ValueError):
        gen.load_holds(hold_file)


def test_two_gene_deposit_needs_a_human(tmp_path):
    """Ruling 4: an intact CDS of the other idiomorph's gene in the window keeps the deposit NEEDS_HUMAN."""
    cache = _copy_cache(tmp_path)
    _edit_genbank(cache, "ORIGIN", "     CDS             2600..2900\n"
                                   "                     /gene=\"MAT1-2-1\"\n"
                                   "                     /product=\"MAT1-2-1\"\n"
                                   "                     /protein_id=\"XXX00001.1\"\n"
                                   "ORIGIN")
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.NEEDS_HUMAN
    assert any(r.startswith("other_idiomorph:") and "XXX00001.1" in r for r in res.reasons("needs_human"))


@pytest.mark.parametrize("location, check", [("join(<1315..1558,1611..2467)", "start:"),
                                             ("join(1315..1558,1611..>2467)", "stop:")])
def test_partial_deposits_need_a_human(tmp_path, location, check):
    """Ruling 10."""
    cache = _copy_cache(tmp_path)
    _edit_genbank(cache, "CDS             join(1315..1558,1611..2467)", f"CDS             {location}")
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.NEEDS_HUMAN
    assert any(r.startswith(check) and "partial" in r for r in res.reasons("needs_human"))


def test_gc_ag_intron_needs_a_human(tmp_path):
    """Ruling 7: a GC-AG intron is a doubtful structure."""
    cache = _copy_cache(tmp_path)
    gbk = _file_starting(cache, "LOCUS       DQ341313")
    lines = gbk.read_text().split("\n")
    i = next(n for n, line in enumerate(lines) if line.strip().startswith("1501 "))
    fields = lines[i].split()
    assert fields[6][8:10] == "gt"
    fields[6] = fields[6][:9] + "c" + fields[6][10:]
    lines[i] = lines[i][:10] + " ".join(fields[1:])
    gbk.write_text("\n".join(lines))
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.NEEDS_HUMAN
    assert any("GC..AG" in r and "doubtful" in r for r in res.reasons("needs_human"))


def test_intron_status_at_the_35_bp_threshold():
    """Ruling 7 at the ruling's own threshold: 34 bp is held, 35 bp passes."""
    assert gen.SHORT_INTRON_BP == 35
    assert gen.intron_status("GT", "AG", 34)[0] == "needs_human"
    assert "shorter than 35 bp" in gen.intron_status("GT", "AG", 34)[1]
    assert gen.intron_status("GT", "AG", 35) == ("pass", "")
    assert gen.intron_status("GC", "AG", 60)[0] == "needs_human"


def test_short_intron_needs_a_human(tmp_path, monkeypatch):
    """Ruling 7 end to end: the deposit's intron is 52 bp; with the threshold above it, it is held."""
    monkeypatch.setattr(gen, "SHORT_INTRON_BP", 60)
    res = _run(_CACHE, _empty_db(tmp_path))
    assert res.verdict == gen.NEEDS_HUMAN
    assert any("shorter than 60 bp" in r for r in res.reasons("needs_human"))


@pytest.mark.parametrize("gene_start", [315, 314])  # 5' extension 1000 bp (caveat) / 1001 bp (held)
def test_wide_gene_feature_threshold(tmp_path, gene_start):
    """Ruling 7: a gene feature reaching > WIDE_GENE_BP (1000) beyond the CDS on a side is held."""
    cache = _copy_cache(tmp_path)
    _edit_genbank(cache, "gene            <1315..>2467", f"gene            {gene_start}..2467")
    res = _run(cache, _empty_db(tmp_path))
    five = 1315 - gene_start
    if five > gen.WIDE_GENE_BP:
        assert res.verdict == gen.NEEDS_HUMAN
        assert any(r.startswith("gene_span:") and "doubtful" in r for r in res.reasons("needs_human"))
    else:
        assert res.verdict == gen.CAVEAT
        assert any(r.startswith("gene_span:") for r in res.reasons("caveat"))
        # ruling 9: the record keeps the CDS span; the gene-feature span is in the note
        assert (res.record["genes"][0]["start"], res.record["genes"][0]["end"]) == (1315, 2467)
        assert f"gene feature is {gene_start}..2467" in res.record["locus"]["core"]["definition_note"]


def test_species_must_be_a_binomial(tmp_path):
    """Ruling 13: organism.species is the species-rank name, and sp./cf. names need a human."""
    for good in ["Ascochyta rabiei", "Melampsora larici-populina"]:
        assert gen.species_binomial(good)
    for bad in ["Xanthoria sp. 1 TBL-2021", "Phialocephala cf. fortinii CSP14", "Beauveria bassiana ARSEF 2860",
                "Mucoraceae sp. AI-2023a", "", None]:
        assert not gen.species_binomial(bad)
    cache = _copy_cache(tmp_path)
    tax = _file_starting(cache, "<?xml")
    text = tax.read_text()
    assert text.count("<ScientificName>Ascochyta rabiei</ScientificName>") == 1
    tax.write_text(text.replace("<ScientificName>Ascochyta rabiei</ScientificName>",
                                "<ScientificName>Ascochyta cf. rabiei</ScientificName>"))
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.NEEDS_HUMAN
    assert any(r.startswith("species:") for r in res.reasons("needs_human"))


def test_disagreeing_curated_record_is_a_conflict(tmp_path):
    """Ruling 8: same protein, different exon model -> NEEDS_HUMAN and a row in the conflicts TSV."""
    repo = tmp_path
    db = _empty_db(repo)
    rec = _run(_CACHE, db).record
    rec["record_id"] = "5454_other-strain_MAT_MAT1-1"
    rec["genes"][0]["exons"] = [{"start": 1300, "end": 1558}, {"start": 1611, "end": 2467}]
    rec["genes"][0]["start"] = 1300
    _write_record(db, rec)
    res = _run(_CACHE, db)
    assert res.verdict == gen.NEEDS_HUMAN
    assert any(r.startswith("conflict:") for r in res.reasons("needs_human"))
    (c,) = res.facts["conflicts"]
    assert c["existing_record_id"] == "5454_other-strain_MAT_MAT1-1"
    assert c["existing_exons"] == "1300-1558;1611-2467" and c["deposit_exons"] == "1315-1558;1611-2467"
    assert c["existing_aa"] == 371 and c["deposit_aa"] == 366 and c["deposit_ncbi_aa"] == 366
    assert "start differs" in c["difference"] and "protein length 371 vs 366" in c["difference"]
    report = tmp_path / "r"
    gen.main(["--repo", str(repo), "--protein", "ABC70320.1", "--date", "2026-10-22", "--cache-dir", str(_CACHE),
              "--offline", "--out", str(tmp_path / "out"), "--report", str(report),
              "--conflicts", str(tmp_path / "conflicts.tsv")])
    rows = (tmp_path / "conflicts.tsv").read_text().strip().split("\n")
    assert len(rows) == 2 and rows[1].startswith("ABC70320.1\tMAT1-1-1\tdb\tsame protein ABC70320.1")
    assert not (tmp_path / "out").exists()  # nothing written


def test_same_strain_identical_model_is_a_duplicate_not_a_conflict(tmp_path):
    db = _empty_db(tmp_path)
    rec = _run(_CACHE, db).record
    _write_record(db, rec)
    res = _run(_CACHE, db)
    assert res.verdict == gen.NEEDS_HUMAN
    assert res.facts["conflicts"] == []
    assert any(r.startswith("duplicate:") and "same model" in r for r in res.reasons("needs_human"))
    assert any(r.startswith("record_id:") for r in res.reasons("needs_human"))  # flagged, not suffixed


# --------------------------------------------------------------------------- fixes after review (2026-10-22)


def _edit_file(cache_dir: Path, prefix: str, old: str, new: str, count: int | None = 1) -> None:
    p = _file_starting(cache_dir, prefix)
    text = p.read_text()
    if count is not None:
        assert text.count(old) == count, (prefix, old)
    p.write_text(text.replace(old, new))


def _as_helotiales(cache: Path) -> None:
    _edit_file(cache, "<?xml", "<ScientificName>Pleosporales</ScientificName>",
               "<ScientificName>Helotiales</ScientificName>")


def _as_mat1_2_1(cache: Path) -> None:
    _edit_file(cache, "LOCUS       DQ341313", "MAT1-1-1", "MAT1-2-1", count=None)
    _edit_file(cache, "LOCUS       ABC70320", "MAT1-1-1", "MAT1-2-1", count=None)


def _without_strain(cache: Path) -> None:
    _edit_file(cache, "LOCUS       DQ341313", '                     /strain="ATCC 76502"\n', "")
    _edit_file(cache, "LOCUS       ABC70320", '                     /strain="ATCC 76502"\n', "")


def _model_record(rec: dict, record_id: str, **gene_changes) -> dict:
    rec = yaml.safe_load(yaml.safe_dump(rec))
    rec["record_id"] = record_id
    rec["genes"][0].update(gene_changes)
    return rec


def test_helotiales_mat1_1_1_is_held_end_to_end(tmp_path):
    """Ruling 1 through build(): the order comes from the taxonomy record."""
    cache = _copy_cache(tmp_path)
    _as_helotiales(cache)
    res = _run(cache, _empty_db(tmp_path), holds=gen.load_holds(None))
    assert res.verdict == gen.NEEDS_HUMAN
    assert res.facts["order"] == "Helotiales"
    assert res.reasons("needs_human") == ["hold: " + gen.HELOTIALES_HOLD_REASON]


def test_helotiales_mat1_2_1_is_not_held(tmp_path):
    cache = _copy_cache(tmp_path)
    _as_helotiales(cache)
    _as_mat1_2_1(cache)
    res = _run(cache, _empty_db(tmp_path), holds=gen.load_holds(None))
    assert res.facts["gene"] == "MAT1-2-1" and res.facts["order"] == "Helotiales"
    assert not any(r.startswith("hold:") for r in res.reasons("needs_human"))
    assert res.verdict == gen.CONFIRMED


def test_hold_file_rows_that_can_match_nothing_warn(tmp_path):
    hold_file = tmp_path / "holds.tsv"
    hold_file.write_text("protein\torder\tgene\treason\n\tErysiphales\tMAT1-1-1\tx\n\tHelotiales\tmat1-1-1\ty\n"
                         "\tHelotiales\tMAT1-1-2\tz\n")
    with pytest.warns(UserWarning) as caught:
        holds = gen.load_holds(hold_file, builtin=False, known_orders={"Helotiales", "Pleosporales"})
    messages = " | ".join(str(w.message) for w in caught)
    assert "order 'Erysiphales'" in messages
    assert "'mat1-1-1' read as 'MAT1-1-1'" in messages
    assert "'MAT1-1-2' is not MAT1-1-1 or MAT1-2-1" in messages
    # lowercase gene is normalised, so it now matches
    assert [h["reason"] for h in gen.matching_holds(holds, "P1.1", "Helotiales", "MAT1-1-1")] == ["y"]


@pytest.mark.parametrize("qualifiers", [
    ['/gene="MAT1-2"', '/product="mating type protein MAT1-2 HMG box"'],
    ['/gene="MAT1-2-1 L"', '/product="putative mating-type 1-2-1 L protein"'],
    ['/product="putative mating-type 1-2-1 R protein"'],
    ['/gene="MAT-2"'],
    ['/gene="mat1-2"'],
])
def test_two_gene_deposit_with_a_loosely_named_other_gene_needs_a_human(tmp_path, qualifiers):
    """Ruling 4 for spellings that are not a roster name (WRK65006.1 / WRK65020.1 use 'MAT1-1-1 L/R')."""
    cache = _copy_cache(tmp_path)
    q = "".join(f"                     {x}\n" for x in qualifiers)
    _edit_genbank(cache, "ORIGIN", f"     CDS             2600..2900\n{q}"
                                   "                     /protein_id=\"XXX00002.1\"\nORIGIN")
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.NEEDS_HUMAN
    assert any(r.startswith("other_idiomorph:") and "XXX00002.1" in r and "MAT1-2" in r
               for r in res.reasons("needs_human"))


def test_a_loosely_named_other_idiomorph_remnant_is_a_caveat_in_the_notes(tmp_path):
    """MT819951.1 layout: a /pseudogene CDS with only /note="dMAT1-x-1"."""
    cache = _copy_cache(tmp_path)
    _edit_genbank(cache, "ORIGIN", "     CDS             2600..2900\n"
                                   "                     /note=\"dMAT1-2-1\"\n"
                                   "                     /pseudogene=\"unprocessed\"\nORIGIN")
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.CAVEAT
    (c,) = [r for r in res.reasons("caveat") if r.startswith("other_idiomorph:")]
    assert "dMAT1-2-1" in c and "pseudo/truncated" in c
    assert "dMAT1-2-1" in res.record["curation"]["notes"]


@pytest.mark.parametrize("gene", ["MAT1-1-2", "MAT1-1-3", "MAT1-1-1 R"])
def test_same_idiomorph_genes_are_not_two_gene_deposits(tmp_path, gene):
    cache = _copy_cache(tmp_path)
    _edit_genbank(cache, "ORIGIN", f"     CDS             2600..2900\n                     /gene=\"{gene}\"\n"
                                   "                     /protein_id=\"XXX00003.1\"\nORIGIN")
    res = _run(cache, _empty_db(tmp_path))
    assert not any(r.startswith("other_idiomorph:") for r in res.reasons("needs_human"))


def test_loose_idiomorph_patterns():
    hit = gen.loose_idiomorphs
    assert hit({"gene": ["MAT1-1-1 L"]}) == {"1"}
    assert hit({"product": ["putative mating-type 1-2-1 protein"]}) == {"2"}
    assert hit({"note": ["dMAT1-1-1"]}) == {"1"}
    assert hit({"gene": ["MAT-2"]}) == {"2"}
    assert hit({"product": ["HMG box protein"]}) == set()
    assert hit({"product": ["hypothetical protein 3"]}) == set()


def test_uncited_deposit_gets_tier_3(tmp_path):
    """Ruling 2, uncited branch."""
    cache = _copy_cache(tmp_path)
    _edit_file(cache, "LOCUS       ABC70320", "   PUBMED   12781674\n", "", count=None)
    res = _run(cache, _empty_db(tmp_path))
    assert res.facts["pmids"] == []
    assert {res.record["evidence"][k]["tier"] for k in res.record["evidence"]} == {3}
    assert "cited2-uncited3 (tier 3)" in res.record["curation"]["notes"]


def test_possible_prediction_is_stated_in_the_notes(tmp_path):
    """Ruling 6, positive sentence."""
    cache = _copy_cache(tmp_path)
    _edit_genbank(cache, '                     /note="mating type protein"\n',
                  '                     /note="mating type protein"\n'
                  '                     /inference="ab initio prediction:AUGUSTUS:3.3"\n')
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.CAVEAT
    notes = res.record["curation"]["notes"]
    assert "ANNOTATION SOURCE: possible gene prediction or third-party annotation" in notes
    assert "AUGUSTUS" in notes


def test_same_strain_on_another_nucleotide_with_other_exon_lengths_is_a_conflict(tmp_path):
    """Ruling 8, taxid + strain criterion and the exon-length path (different accession)."""
    db = _empty_db(tmp_path)
    rec = _run(_CACHE, db).record
    other = _model_record(rec, "5454_atcc-76502_MAT_MAT1-1", protein_accession="ncbi_protein:ZZZ00001.1",
                          exons=[{"start": 100, "end": 400}, {"start": 460, "end": 1316}])
    other["locus"]["core"]["segments"][0]["sequence_source"].update(accession="ZZ999999.1", seq_region="ZZ999999.1")
    _write_record(db, other)
    res = _run(_CACHE, db)
    assert res.verdict == gen.NEEDS_HUMAN
    (c,) = res.facts["conflicts"]
    assert c["source"] == "db" and c["match"] == "same taxid 5454 and strain atcc-76502"
    assert "exon lengths [301, 857] vs [244, 857]" in c["difference"]


def test_same_taxid_other_strain_other_nucleotide_is_not_a_conflict(tmp_path):
    """Negative control: a matcher on taxid alone would fail this."""
    db = _empty_db(tmp_path)
    rec = _run(_CACHE, db).record
    other = _model_record(rec, "5454_cbs-1_MAT_MAT1-1", protein_accession="ncbi_protein:ZZZ00001.1",
                          exons=[{"start": 100, "end": 400}, {"start": 460, "end": 1316}])
    other["organism"]["strain"]["name"] = "CBS 1"
    other["locus"]["core"]["segments"][0]["sequence_source"].update(accession="ZZ999999.1", seq_region="ZZ999999.1")
    _write_record(db, other)
    res = _run(_CACHE, db)
    assert res.facts["conflicts"] == []
    assert res.verdict == gen.CAVEAT
    assert any(r.startswith("same_species:") and "db/ already holds" in r for r in res.reasons("caveat"))


def test_same_strain_under_a_strain_level_taxid_is_matched(tmp_path):
    """Ruling 8: curated record filed under a strain-level taxid, deposit under the species taxid."""
    db = _empty_db(tmp_path)
    rec = _run(_CACHE, db).record
    other = _model_record(rec, "999999_atcc-76502_MAT_MAT1-1", protein_accession="ncbi_protein:ZZZ00001.1",
                          exons=[{"start": 100, "end": 400}, {"start": 460, "end": 1316}])
    other["taxonomy"]["taxid"] = 999999
    other["locus"]["core"]["segments"][0]["sequence_source"].update(accession="ZZ999999.1", seq_region="ZZ999999.1")
    _write_record(db, other)
    res = _run(_CACHE, db)
    (c,) = res.facts["conflicts"]
    assert c["match"] == "same species Ascochyta rabiei (taxid 999999 vs 5454) and strain atcc-76502"


def test_assembly_segment_seq_region_is_indexed(tmp_path):
    """Ruling 8: an assembly-typed segment keeps the contig in seq_region."""
    db = _empty_db(tmp_path)
    rec = _run(_CACHE, db).record
    other = _model_record(rec, "5454_cbs-2_MAT_MAT1-1", protein_accession="ncbi_protein:ZZZ00001.1",
                          exons=[{"start": 1300, "end": 1558}, {"start": 1611, "end": 2467}], start=1300)
    other["locus"]["core"]["segments"][0]["sequence_source"] = {
        "type": "assembly", "accession": "GCA_000000001.1", "seq_region": "DQ341313.1"}
    _write_record(db, other)
    res = _run(_CACHE, db)
    (c,) = res.facts["conflicts"]
    assert c["match"] == "overlapping region of DQ341313.1" and "start differs" in c["difference"]


def test_existing_gene_without_exons_and_equal_span_is_not_a_conflict(tmp_path):
    """Ruling 8: a span-only curated gene whose span equals the CDS is the same gene (duplicate)."""
    db = _empty_db(tmp_path)
    rec = _run(_CACHE, db).record
    other = _model_record(rec, "5454_atcc-76502_MAT_MAT1-1", exons=[], strand=None)
    del other["genes"][0]["exons"]
    _write_record(db, other)
    res = _run(_CACHE, db)
    assert res.facts["conflicts"] == []
    assert any(r.startswith("duplicate:") and "lacks exons; its span equals" in r for r in res.reasons("needs_human"))


def test_existing_gene_without_exons_on_other_coordinates_is_uncompared(tmp_path):
    db = _empty_db(tmp_path)
    rec = _run(_CACHE, db).record
    other = _model_record(rec, "5454_atcc-76502_MAT_MAT1-1", protein_accession="ncbi_protein:ZZZ00001.1",
                          start=10, end=1162)
    del other["genes"][0]["exons"]
    other["locus"]["core"]["segments"][0]["sequence_source"].update(accession="ZZ999999.1", seq_region="ZZ999999.1")
    _write_record(db, other)
    res = _run(_CACHE, db)
    assert res.facts["conflicts"] == []
    assert any(r.startswith("existing_no_exons:") for r in res.reasons("needs_human"))


def test_existing_identical_model_with_null_strand_is_not_a_conflict(tmp_path):
    db = _empty_db(tmp_path)
    rec = _run(_CACHE, db).record
    _write_record(db, _model_record(rec, "5454_atcc-76502_MAT_MAT1-1", strand=None))
    res = _run(_CACHE, db)
    assert res.facts["conflicts"] == []
    assert any(r.startswith("duplicate:") and "same model" in r for r in res.reasons("needs_human"))


def test_in_batch_candidates_are_not_called_curated(tmp_path):
    """The same deposit twice in one run: the second is a batch duplicate, never 'curated'."""
    repo = tmp_path
    _empty_db(repo)
    report = tmp_path / "r"
    gen.main(["--repo", str(repo), "--protein", "ABC70320.1", "--protein", "ABC70320.1", "--date", "2026-10-22",
              "--cache-dir", str(_CACHE), "--offline", "--report", str(report)])
    lines = Path(f"{report}.tsv").read_text().strip().split("\n")[1:]
    first, second = (line.split("\t") for line in lines)
    assert first[1] == gen.CONFIRMED and second[1] == gen.NEEDS_HUMAN
    assert "batch_duplicate:" in second[12] and "in-batch candidate (not written)" in second[12]
    assert "curated" not in second[12]


def test_in_batch_candidates_never_change_unknown_slugs(tmp_path):
    cache = _copy_cache(tmp_path)
    _without_strain(cache)
    db = _empty_db(tmp_path)
    batch = gen.new_index("batch")
    gen.add_record_to_index(batch, {"record_id": "5454_unknown-1_MAT_MAT1-2", "taxonomy": {"taxid": 5454}})
    res = _run(cache, db, batch_index=batch)
    assert res.facts["record_id"] == "5454_unknown-1_MAT_MAT1-1"
    assert res.verdict != gen.NEEDS_HUMAN
    # a curated db/ record does advance the counter (the spec rule)
    rec = yaml.safe_load(yaml.safe_dump(res.record))
    rec["record_id"] = "5454_unknown-1_MAT_MAT1-2"
    rec["genes"][0]["protein_accession"] = "ncbi_protein:ZZZ00009.1"
    rec["locus"]["core"]["segments"][0]["sequence_source"].update(accession="ZZ999998.1", seq_region="ZZ999998.1")
    _write_record(db, rec)
    assert _run(cache, db).facts["record_id"] == "5454_unknown-2_MAT_MAT1-1"


def test_other_idiomorph_accessory_gene_is_held_with_its_own_reason(tmp_path):
    cache = _copy_cache(tmp_path)
    _edit_genbank(cache, "ORIGIN", "     CDS             2600..2900\n                     /product=\"MAT1-2-4\"\n"
                                   "                     /protein_id=\"XXX00004.1\"\nORIGIN")
    res = _run(cache, _empty_db(tmp_path))
    assert res.verdict == gen.NEEDS_HUMAN
    assert any("accessory gene of the other idiomorph" in r and "XXX00004.1" in r for r in res.reasons("needs_human"))


def test_a_locus_note_on_a_neighbouring_orf_is_not_a_gene_name(tmp_path):
    """AKE48510.1 layout: a hypothetical protein with /note="MAT1-2 mating-type locus"."""
    cache = _copy_cache(tmp_path)
    _edit_genbank(cache, "ORIGIN", "     CDS             2600..2900\n"
                                   "                     /product=\"hypothetical protein\"\n"
                                   "                     /note=\"MAT1-2 mating-type locus\"\n"
                                   "                     /protein_id=\"XXX00005.1\"\nORIGIN")
    res = _run(cache, _empty_db(tmp_path))
    assert not any(r.startswith("other_idiomorph:") for r in res.reasons("needs_human") + res.reasons("caveat"))
