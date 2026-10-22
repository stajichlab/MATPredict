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
    # the schema does not accept `unknown` yet (its own PR); that is the record's only schema error
    errors = validate_metadata(res.record)
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


def test_short_intron_needs_a_human(tmp_path, monkeypatch):
    """Ruling 7: the deposit's intron is 52 bp; with the threshold above it, it is held."""
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
    assert len(rows) == 2 and rows[1].startswith("ABC70320.1\tMAT1-1-1\tsame protein ABC70320.1")
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
