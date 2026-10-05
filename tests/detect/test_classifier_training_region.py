"""A record gene's `classifier_training_region` trims what the classifier trains on.

Curator ruling 2026-10-04: the Mycotypha africana sexM is a 531-aa ORF whose
C-terminal HMG-like repeats should not enter the HMM; the record keeps the full
protein (and search uses it), training uses the HMG-region segment only.
"""
from pathlib import Path

import yaml

from MATPredict.detect import classifier_build as cb


def _write_record(root: Path, record_id: str, region=None):
    d = root / "Mucoromycota" / "Mucorales" / record_id
    d.mkdir(parents=True)
    gene = {"gene_index": 0, "name": "sexM", "role": "core_MAT", "present": True}
    if region:
        gene["classifier_training_region"] = {"aa_start": region[0], "aa_end": region[1], "reason": "test"}
    (d / "metadata.yaml").write_text(yaml.safe_dump({"record_id": record_id, "genes": [gene]}))


def test_regions_are_read_from_accepted_records_only(tmp_path):
    _write_record(tmp_path, "R1", (26, 218))
    _write_record(tmp_path, "R2")
    cand = tmp_path / "candidates" / "Mucoromycota" / "R3"
    cand.mkdir(parents=True)
    (cand / "metadata.yaml").write_text(yaml.safe_dump({"record_id": "R3", "genes": [
        {"gene_index": 0, "name": "sexM", "classifier_training_region": {"aa_start": 1, "aa_end": 5, "reason": "x"}}]}))
    assert cb._training_regions(tmp_path) == {("R1", 0): (26, 218)}


def test_training_set_trims_to_the_region(tmp_path, monkeypatch):
    _write_record(tmp_path, "R1", (3, 62))
    full = "M" + "A" * 99  # 100 aa
    other = "M" + "C" * 99

    def fake_ref(db_root, out, family_keys=None):
        Path(out).write_text(f">R1|gene0|sexM\n{full}\n>R2|gene0|sexM\n{other}\n")
        return Path(out)

    class Fam:
        key = "Mucoromycota:MAT"

    monkeypatch.setattr(cb, "build_reference_fasta", fake_ref)
    monkeypatch.setattr(cb, "classifier_genes", lambda fam: {"sexM": "Minus"})
    rows = {r["id"]: r["sequence"] for r in cb.training_set(tmp_path, Fam(), tmp_path / "cls")}
    assert rows["REF|R1|gene0|sexM"] == full[2:62]
    assert rows["REF|R2|gene0|sexM"] == other
