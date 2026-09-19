from __future__ import annotations

import csv

import pytest
import yaml

from MATPredict.db.gff_export import write_genbank
from MATPredict.db.synteny import ClinkerError, SyntenyRecordError, draw_synteny

# Two real curated-record shapes sharing gene names across records -- the real,
# common case for this project (e.g. COX13/APN2/MAT1-1-1 recur across every
# curated record of the same locus family; see db/Ascomycota/Onygenales's real
# Coccidioides MAT1-1/MAT1-2 records) -- plus a third record with a distinct
# gene, so the fixture set exercises both a shared-name role/colour mapping and
# a per-record-unique one.
RECORD_A = {
    "record_id": "111_a_MAT_MAT1-1",
    "locus": {"core": {"segments": [
        {"segment_index": 0, "start": 1, "end": 300, "sequence_source": {"seq_region": "scaffold_1"}},
    ]}},
    "genes": [
        {"gene_index": 0, "name": "COX13", "role": "flanking_conserved", "present": True,
         "segment_index": 0, "start": 10, "end": 100, "strand": "+"},
        {"gene_index": 1, "name": "MAT1-1-1", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 110, "end": 200, "strand": "+"},
    ],
}

RECORD_B = {
    "record_id": "222_b_MAT_MAT1-2",
    "locus": {"core": {"segments": [
        {"segment_index": 0, "start": 1, "end": 300, "sequence_source": {"seq_region": "scaffold_2"}},
    ]}},
    "genes": [
        {"gene_index": 0, "name": "COX13", "role": "flanking_conserved", "present": True,
         "segment_index": 0, "start": 10, "end": 100, "strand": "+"},
        {"gene_index": 1, "name": "MAT1-2-1", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 110, "end": 200, "strand": "+"},
        {"gene_index": 2, "name": "SLA2like", "role": "flanking_variable", "present": True,
         "segment_index": 0, "start": 210, "end": 280, "strand": "+",
         "gene_class": "flanking_extended", "present_in_idiomorphs": ["MAT1-2"]},
        # A curated-but-absent gene: never a real row in the generated CSVs.
        {"gene_index": 3, "name": "GhostGene", "role": "core_MAT", "present": False,
         "segment_index": 0, "start": 290, "end": 295, "strand": "+"},
    ],
}


def _write_record(tmp_path, record):
    """Write one accepted record's real locus.gbk + metadata.yaml under a
    db_root, matching this project's db/<Phylum>/<Order>/<record_id>/ layout."""
    record_dir = tmp_path / "Ascomycota" / "Onygenales" / record["record_id"]
    record_dir.mkdir(parents=True)
    sequences = {gene["gene_index"]: "M" * 20 for gene in record["genes"] if gene.get("present", True)}
    write_genbank(record, sequences, out_path=record_dir / "locus.gbk")
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(record, sort_keys=False))
    return record_dir


class _Result:
    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def test_draw_synteny_invokes_clinker_with_resolved_gbk_paths_and_csvs(tmp_path):
    db_root = tmp_path / "db"
    dir_a = _write_record(db_root, RECORD_A)
    dir_b = _write_record(db_root, RECORD_B)
    out_path = tmp_path / "out" / "synteny.html"

    commands = []

    def fake_runner(cmd, **kwargs):
        commands.append(cmd)
        return _Result()

    result = draw_synteny(
        [RECORD_A["record_id"], RECORD_B["record_id"]], db_root, out_path, runner=fake_runner
    )

    assert result == out_path
    assert len(commands) == 1
    cmd = commands[0]
    assert cmd[0] == "clinker"
    assert str(dir_a / "locus.gbk") in cmd
    assert str(dir_b / "locus.gbk") in cmd
    assert "-gf" in cmd
    assert "-cm" in cmd
    assert "-p" in cmd
    assert cmd[cmd.index("-p") + 1] == str(out_path)
    assert "-f" in cmd


def test_draw_synteny_generates_gene_functions_and_colour_map_csvs(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, RECORD_A)
    _write_record(db_root, RECORD_B)
    out_path = tmp_path / "out" / "synteny.html"

    def fake_runner(cmd, **kwargs):
        return _Result()

    draw_synteny([RECORD_A["record_id"], RECORD_B["record_id"]], db_root, out_path, runner=fake_runner)

    gene_functions_path = out_path.parent / "synteny.gene_functions.csv"
    colour_map_path = out_path.parent / "synteny.colour_map.csv"
    assert gene_functions_path.exists()
    assert colour_map_path.exists()

    with gene_functions_path.open() as fh:
        rows = list(csv.reader(fh))
    # gene_id -> role, exactly the real curated genes across both records, in the
    # gene qualifier clinker itself will key on ("gene" -> gene["name"]); no header
    # row (parse_gene_functions has no header-skip, and would treat one as data);
    # the not-present GhostGene must never appear.
    assert rows == [
        ["COX13", "flanking_conserved"],
        ["MAT1-1-1", "core_MAT"],
        ["COX13", "flanking_conserved"],
        ["MAT1-2-1", "core_MAT"],
        ["SLA2like", "flanking_variable"],
    ]

    with colour_map_path.open() as fh:
        colour_rows = dict(csv.reader(fh))
    # Reuses draw.py's exact ROLE_COLORS palette for the 3 real roles present.
    assert colour_rows == {
        "core_MAT": "#D55E00",
        "flanking_conserved": "#0072B2",
        "flanking_variable": "#009E73",
    }


def test_draw_synteny_requires_at_least_two_records(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, RECORD_A)
    with pytest.raises(ValueError):
        draw_synteny([RECORD_A["record_id"]], db_root, tmp_path / "out.html")


def test_draw_synteny_reports_a_clear_error_for_an_unresolvable_record_id(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, RECORD_A)
    with pytest.raises(SyntenyRecordError, match="no_such_record"):
        draw_synteny([RECORD_A["record_id"], "no_such_record"], db_root, tmp_path / "out.html")


def test_draw_synteny_ignores_candidate_records(tmp_path):
    """A record under db_root/candidates/... is not accepted and must not be
    silently used even if its record_id happens to match."""
    db_root = tmp_path / "db"
    candidate_dir = db_root / "candidates" / "Ascomycota" / RECORD_A["record_id"]
    candidate_dir.mkdir(parents=True)
    sequences = {gene["gene_index"]: "M" * 20 for gene in RECORD_A["genes"]}
    write_genbank(RECORD_A, sequences, out_path=candidate_dir / "locus.gbk")
    (candidate_dir / "metadata.yaml").write_text(yaml.safe_dump(RECORD_A, sort_keys=False))
    _write_record(db_root, RECORD_B)

    with pytest.raises(SyntenyRecordError, match=RECORD_A["record_id"]):
        draw_synteny([RECORD_A["record_id"], RECORD_B["record_id"]], db_root, tmp_path / "out.html")


def test_draw_synteny_raises_clinker_error_on_nonzero_returncode(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, RECORD_A)
    _write_record(db_root, RECORD_B)

    def failing_runner(cmd, **kwargs):
        return _Result(returncode=1, stderr="clinker: no alignments found")

    with pytest.raises(ClinkerError, match="clinker"):
        draw_synteny(
            [RECORD_A["record_id"], RECORD_B["record_id"]], db_root, tmp_path / "out.html",
            runner=failing_runner,
        )
