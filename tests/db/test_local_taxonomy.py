"""Offline NCBI taxonomy table (`MATPredict.db.local_taxonomy`).

Curator ruling, 2026-10-04: embed a slim all-taxa table in the container so a
run needs no NCBI E-utilities call.

Measured on the real NCBI archive taxdmp_2026-10-01.zip: 3,016,750 taxa,
101,411 merged ids, 23.7 MB table; load 1.96 s and 76 MB RSS. Against 8,185
taxonomy documents in the run cache (efetch answers fetched 2026-09/10):
phylum 8,185/8,185 and genetic code 8,185/8,185 identical; lineage 8,179/8,185
(the 6 differ because NCBI re-parented the taxa between the fetch and the
snapshot, e.g. C. gattii VGV 2268382 lost the no-rank node 3111917); 21 taxids
are deleted in the snapshot (delnodes.dmp).
"""
import gzip
import zipfile

import pytest

from MATPredict.db import local_taxonomy as lt


# root(1) - cellular(131567) - Eukaryota(2759) - Fungi(4751, kingdom)
#   - Ascomycota(4890, phylum) - Pachysolen(4917, genus) - P. tannophilus(4918, code 26)
#   - Mucoromycota(1913637, phylum) - Mucor(4830, genus) - M. circinelloides(4837)
NODES = [
    (1, 1, "no rank", 1), (131567, 1, "cellular root", 1), (2759, 131567, "domain", 1),
    (4751, 2759, "kingdom", 1), (4890, 4751, "phylum", 1), (4917, 4890, "genus", 26),
    (4918, 4917, "species", 26), (1913637, 4751, "phylum", 1), (4830, 1913637, "genus", 1),
    (4837, 4830, "species", 1),
]
NAMES = {1: "root", 131567: "cellular organisms", 2759: "Eukaryota", 4751: "Fungi",
         4890: "Ascomycota", 4917: "Pachysolen", 4918: "Pachysolen tannophilus",
         1913637: "Mucoromycota", 4830: "Mucor", 4837: "Mucor circinelloides"}


def _dmp(rows):
    return "".join("\t|\t".join(str(x) for x in r) + "\t|\n" for r in rows)


def _dump_texts():
    nodes = _dmp([(t, p, r, "", 0, 1, c, 1, 1, 1, 0, 0, "") for t, p, r, c in NODES])
    names = _dmp([(t, n, "", "scientific name") for t, n in NAMES.items()]
                 + [(4837, "Mucor racemosus f. circinelloides", "", "synonym")])
    merged = _dmp([(99999, 4837)])
    return {"nodes.dmp": nodes, "names.dmp": names, "merged.dmp": merged}


@pytest.fixture
def table(tmp_path):
    d = tmp_path / "dump"
    d.mkdir()
    for name, text in _dump_texts().items():
        (d / name).write_text(text)
    out = tmp_path / "tax.tsv.zst"
    counts = lt.build_table(d, out, snapshot="2026-10-01")
    assert counts == {"taxa": len(NODES), "merged": 1, "snapshot": "2026-10-01"}
    return out


def test_lookups_match_efetch_semantics(table):
    t = lt.LocalTaxonomy(table)
    assert t.snapshot == "2026-10-01"
    # efetch LineageEx: root-first, no root (1), not the taxid itself
    assert t.lineage(4918) == [131567, 2759, 4751, 4890, 4917]
    assert t.phylum_name(4918) == "Ascomycota"
    assert t.genetic_code(4918) == 26
    assert t.phylum_name(4837) == "Mucoromycota"
    # a phylum's own lineage has no phylum ancestor (efetch gives none either)
    assert t.phylum_name(4890) is None


def test_merged_taxid_answers_for_its_new_taxid(table):
    t = lt.LocalTaxonomy(table)
    assert t.lineage(99999) == t.lineage(4837)
    assert 99999 in t


def test_unknown_taxid_raises(table):
    t = lt.LocalTaxonomy(table)
    assert 123456 not in t
    with pytest.raises(lt.TaxidNotInSnapshot, match="2026-10-01"):
        t.lineage(123456)


@pytest.mark.parametrize("form", ["gz_dir", "zip"])
def test_builds_from_compressed_dumps(tmp_path, form):
    texts = _dump_texts()
    if form == "gz_dir":
        src = tmp_path / "d"
        src.mkdir()
        for name, text in texts.items():
            (src / (name + ".gz")).write_bytes(gzip.compress(text.encode()))
    else:
        src = tmp_path / "taxdmp_2026-10-01.zip"
        with zipfile.ZipFile(src, "w") as zf:
            for name, text in texts.items():
                zf.writestr(name, text)
    out = tmp_path / "t.tsv.zst"
    lt.build_table(src, out, snapshot="x")
    assert lt.LocalTaxonomy(out).genetic_code(4918) == 26


def test_default_resolvers_use_the_table_without_network(table, monkeypatch):
    from MATPredict.db import taxonomy
    monkeypatch.setenv("MATPREDICT_TAXONOMY", str(table))
    monkeypatch.setattr(taxonomy, "_get_default_ncbi_client",
                        lambda: (_ for _ in ()).throw(AssertionError("network used")))
    assert taxonomy.default_lineage_taxids(4918)[-1] == 4917
    assert taxonomy.default_lineage_phylum_name(4918) == "Ascomycota"
    assert taxonomy.default_genetic_code(4918) == 26


def test_offline_unknown_taxid_raises_instead_of_calling_ncbi(table, monkeypatch):
    from MATPredict.db import taxonomy
    monkeypatch.setenv("MATPREDICT_TAXONOMY", str(table))
    monkeypatch.setenv("MATPREDICT_OFFLINE", "1")
    monkeypatch.setattr(taxonomy, "_get_default_ncbi_client",
                        lambda: (_ for _ in ()).throw(AssertionError("network used")))
    with pytest.raises(lt.TaxidNotInSnapshot, match="--phylum"):
        taxonomy.default_lineage_taxids(123456)


def test_online_unknown_taxid_falls_back_to_ncbi(table, monkeypatch):
    from MATPredict.db import taxonomy
    monkeypatch.setenv("MATPREDICT_TAXONOMY", str(table))
    monkeypatch.delenv("MATPREDICT_OFFLINE", raising=False)

    class Client:
        def fetch_taxonomy_lineage(self, taxid):
            return [42]
    monkeypatch.setattr(taxonomy, "_get_default_ncbi_client", lambda: Client())
    assert taxonomy.default_lineage_taxids(123456) == [42]


def test_offline_without_a_table_raises(monkeypatch):
    from MATPredict.db import taxonomy
    monkeypatch.delenv("MATPREDICT_TAXONOMY", raising=False)
    monkeypatch.setenv("MATPREDICT_OFFLINE", "1")
    with pytest.raises(lt.TaxidNotInSnapshot, match="no local taxonomy table"):
        taxonomy.default_genetic_code(4918)


def test_source_label(table, monkeypatch):
    monkeypatch.delenv("MATPREDICT_TAXONOMY", raising=False)
    monkeypatch.delenv("MATPREDICT_OFFLINE", raising=False)
    assert lt.source_label() == "NCBI E-utilities"
    monkeypatch.setenv("MATPREDICT_TAXONOMY", str(table))
    assert lt.source_label().startswith("local NCBI taxonomy snapshot 2026-10-01")
    monkeypatch.setenv("MATPREDICT_OFFLINE", "1")
    assert lt.source_label() == "local NCBI taxonomy snapshot 2026-10-01"
