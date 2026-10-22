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


def _run(cache_dir: Path):
    fetcher = gen.Fetcher(cache_dir, offline=True)
    db_root = _ROOT / "db"
    return gen.build("ABC70320.1", fetcher, db_root, gen.load_db_index(db_root), date="2026-10-22",
                     mating_system="heterothallic", tier_policy="cited2-uncited3")


def _copy_cache(tmp_path: Path) -> Path:
    dst = tmp_path / "cache"
    shutil.copytree(_CACHE, dst)
    return dst


def _file_starting(cache_dir: Path, prefix: str) -> Path:
    return next(p for p in cache_dir.iterdir() if p.read_text().startswith(prefix))


def test_real_deposit_gives_a_valid_confirmed_record():
    res = _run(_CACHE)
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
    # the record passes the same checks `curate-db propose/accept/validate` apply
    order_doc = yaml.safe_load((_ROOT / "db" / "Ascomycota" / "order.yml").read_text())
    assert validate_metadata(rec) == []
    assert validate_idiomorphs(rec, order_doc) == []
    assert validate_gene_vocabulary(rec, order_doc) == []
    assert validate_gene_status(rec) == []


def test_output_is_deterministic():
    a = yaml.safe_dump(_run(_CACHE).record, sort_keys=False)
    b = yaml.safe_dump(_run(_CACHE).record, sort_keys=False)
    assert a == b


def test_a_protein_that_the_exons_do_not_encode_needs_a_human(tmp_path):
    cache = _copy_cache(tmp_path)
    fasta = _file_starting(cache, ">ABC70320.1")
    header, seq = fasta.read_text().split("\n", 1)
    seq = seq.replace("\n", "")
    mutated = seq[:20] + ("A" if seq[20] != "A" else "G") + seq[21:]
    fasta.write_text(header + "\n" + mutated + "\n")
    res = _run(cache)
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
    res = _run(cache)
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
