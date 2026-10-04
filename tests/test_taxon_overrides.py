"""db/taxon_overrides.tsv and its loader (curator ruling 2026-10-01, B8)."""
from __future__ import annotations

from pathlib import Path

import pytest

from MATPredict.detect.taxon_overrides import COLUMNS, load_taxon_overrides

DB = Path(__file__).resolve().parents[1] / "db"
HEADER = "\t".join(COLUMNS)


def test_shipped_file_loads_the_ruled_strains():
    ov = load_taxon_overrides(DB)
    assert set(ov) == {"Backusella_ctenidia_NRRL_6239", "Rhizopus_arrhizus_NRRL_1470",
                       "Thamnostylum_repens_Tieghem_Upadhyay_NRRL_6240", "Rhizopus_microsporus_NRRL_A-17693",
                       "Absidia_sp._NRRL_3163",
                       # B12 ITS ruling 2026-10-02: tree genus and rDNA agree (confirmed),
                       # except Circinella_muscae_NRRL_1360 (rDNA and markers disagree).
                       "Phycomyces_blakesleeanus_NRRL_1555", "Phycomyces_blakesleeanus_NRRL_1556",
                       "Phycomyces_nitens_NRRL_2700", "Pilaira_anomala_NRRL_2527",
                       "Pilaira_anomala_RSA_1997_Plus", "Thamnidium_elegans_NRRL_2467",
                       "Syncephalastrum_racemosum_NRRL_1506", "Syncephalastrum_racemosum_NRRL_1623",
                       "Syncephalastrum_sp._NRRL_1485", "Rhizomucor_pusillus_NRRL_2543",
                       "Cunninghamella_japonica_NRRL_2464", "Circinella_tenella_NRRL_A-23557",
                       "Circinella_muscae_NRRL_1360", "Mucor_sp._NRRL_1454"}
    assert ov["Rhizopus_microsporus_NRRL_A-17693"].likely_identity == "Circinella minor"
    confirmed = {g for g, o in ov.items() if o.status == "confirmed"}
    assert len(confirmed) == 13
    assert ov["Circinella_muscae_NRRL_1360"].status == "unconfirmed"


def test_missing_file_is_empty(tmp_path):
    assert load_taxon_overrides(tmp_path) == {}


def _write(tmp_path, *rows):
    p = tmp_path / "taxon_overrides.tsv"
    p.write_text("# comment\n" + HEADER + "\n" + "\n".join(rows) + "\n")
    return p


ROW = "g1\tA b\tC d\tspecies\tidentical locus\tunconfirmed\texclude\tpath\tcurator 2026-10-01"


def test_valid_row(tmp_path):
    assert load_taxon_overrides(_write(tmp_path, ROW))["g1"].identity_rank == "species"


@pytest.mark.parametrize("row,msg", [
    (ROW.replace("unconfirmed", "maybe"), "status"),
    (ROW.replace("species", "order"), "identity_rank"),
    (ROW.replace("identical locus", ""), "empty basis"),
    ("g1\tA b", "fields"),
])
def test_malformed_rows_raise(tmp_path, row, msg):
    with pytest.raises(ValueError, match=msg):
        load_taxon_overrides(_write(tmp_path, row))


def test_duplicate_raises(tmp_path):
    with pytest.raises(ValueError, match="duplicate"):
        load_taxon_overrides(_write(tmp_path, ROW, ROW))
