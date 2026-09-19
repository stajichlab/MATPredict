"""Multi-locus synteny/comparison diagram across 2+ curated records, via clinker.

This is a curated-record-side (`db/`) feature, complementary to `draw.py`'s
single-locus pyGenomeViz figure: where `draw.py` renders one record's own gene
structure, this module compares several records' `locus.gbk` files
(Task 1's `gff_export.write_genbank` output) side by side, auto-generating
clinker's `--gene_functions`/`--colour_map` CSV inputs directly from each
gene's already-curated `role` schema field, rather than requiring a curator to
hand-build those CSVs.

Real, confirmed clinker CLI shape (verified against the installed
`clinker-py` 0.0.32 this session, both via `pixi run clinker --help` and by
reading the installed package source directly, per Task 3's own verification
discipline)::

    clinker <genbank files...> -gf/--gene_functions <csv> -cm/--colour_map <csv> \
        -p/--plot <html_path> -f/--force

Column semantics (`clinker/main.py`'s real `parse_gene_functions`/
`parse_colour_map`, and `clinker/align.py`'s real `Globaligner.build_gene_groups`/
`get_gene_uid`):

- `--gene_functions` is a HEADERLESS 2-column CSV, `gene_id,function_label`.
  `parse_gene_functions` reads it with the plain stdlib `csv` module and
  inverts it into `{function_label: [gene_id, ...]}` -- so a header row would
  itself be read as a bogus `("Gene", "Function")` data row and must not be
  written.
- `--colour_map` is a SEPARATE headerless 2-column CSV, `function_label,hex`.
  Its first column must use the SAME label strings `--gene_functions`' second
  column uses -- it keys on the function label, not on any gene identifier.
- `get_gene_uid` resolves each `gene_id` string in `--gene_functions` against
  every loaded gene's own `label` property, which clinker derives via
  `find_qualifier` over each GenBank feature's qualifiers, in this fixed
  precedence: `protein_id`, `locus_tag`, `id`, `ID`, `gene`, `label`, `name`.
  `gff_export.write_genbank`'s `gene`/`CDS` features carry none of
  `protein_id`/`locus_tag`/`id`/`ID`/`label` -- only `gene` -- so clinker's
  real label for every one of this project's genes is exactly the gene's
  `name` schema field. `--gene_functions`' first column therefore must be
  each gene's `name`, unmodified.

Known, confirmed clinker limitation this module cannot work around (verified
by reading `Globaligner.get_gene_uid`/`get_gene_uids`, not guessed): `_genes`
is keyed by an internal per-feature uid, but `get_gene_uid(label)` does a
linear search over `_genes.items()` and returns the FIRST uid whose stored
label matches. When two DIFFERENT genes across two DIFFERENT input records
share the same curated `name` (a common, expected case for this project --
e.g. `MAT1-1-1`/`COX13`/`APN2` recur across every curated record of the same
locus family), only the first-encountered occurrence of that name is ever
added to its function group; later occurrences of the same name are not
colored/grouped by function even though `--gene_functions` lists them too.
This is a real limitation of clinker's own CLI, not introduced by this
module's CSV generation -- `--gene_functions`' gene-id column has no way to
distinguish two features that clinker itself will resolve to the identical
label. Confirmed live in this task's own end-to-end trial (see
`docs/superpowers/sdd/2026-09-18-mat-locus-visualization/task-5-report.md`).

Function-label design choice: the function label is the gene's `role` value
(`core_MAT` / `flanking_conserved` / `flanking_variable`) alone, not
`role:gene_class`. Folding `gene_class` into the label would fragment genes
that share a role but differ in `gene_class` into separate colour groups,
losing the at-a-glance "which genes are core vs flanking" comparison this
diagram exists for, and would require the colour map to carry one row per
distinct `role`+`gene_class` combination rather than reusing `draw.py`'s
fixed 3-color `ROLE_COLORS` palette (kept identical here for visual
consistency between the two diagram types, per this project's own visual
convention).
"""
from __future__ import annotations

import csv
import subprocess
from pathlib import Path
from typing import Callable

import yaml

from MATPredict.db.draw import ROLE_COLORS

_FALLBACK_COLOR = "#999999"


class SyntenyRecordError(ValueError):
    """A requested record id did not resolve to exactly one curated locus.gbk."""


class ClinkerError(RuntimeError):
    """The real `clinker` CLI exited non-zero."""


def _run_checked(runner: Callable, cmd: list[str]):
    """Run cmd and raise ClinkerError if it exits non-zero.

    Mirrors `detect.search._run_checked`'s convention: without this check, a
    missing binary or a malformed input would produce a silent, wrong (or
    default-colored) diagram rather than a clear failure.
    """
    result = runner(cmd, capture_output=True, text=True)
    returncode = getattr(result, "returncode", 0)
    if returncode != 0:
        stderr = (getattr(result, "stderr", "") or "").strip()
        raise ClinkerError(
            f"`{' '.join(cmd)}` exited {returncode}" + (f": {stderr}" if stderr else "")
        )
    return result


def _resolve_record(record_id: str, db_root: Path) -> tuple[Path, Path]:
    """Resolve a record id to its real (locus.gbk, metadata.yaml) pair under db_root.

    Matches this project's existing `db/<Phylum>/<Order>/<record_id>/` layout
    (the same 3-level shape `db/cli.py`'s `find_records_missing_proteins_faa`
    walks), explicitly excluding `db_root/candidates/...` -- only an accepted
    record has a `locus.gbk` clinker can actually use.
    """
    matches = sorted(
        p for p in db_root.glob(f"*/*/{record_id}/locus.gbk")
        if p.relative_to(db_root).parts[0] != "candidates"
    )
    if not matches:
        raise SyntenyRecordError(
            f"no accepted record {record_id!r} with a locus.gbk found under {db_root} "
            "(the record id may not exist, may still be a candidate, or may not have "
            "had `matpredict curate-db build-gff` run for it yet)"
        )
    if len(matches) > 1:
        raise SyntenyRecordError(
            f"record id {record_id!r} is ambiguous under {db_root}: matches {matches}"
        )
    gbk_path = matches[0]
    metadata_path = gbk_path.parent / "metadata.yaml"
    if not metadata_path.exists():
        raise SyntenyRecordError(f"{gbk_path} has no sibling metadata.yaml at {metadata_path}")
    return gbk_path, metadata_path


def _gene_function_rows(metadata_path: Path) -> list[tuple[str, str]]:
    """(gene name, role) for every present gene in one record's real metadata.yaml.

    Never fabricates a role: a present gene with no `role` field is skipped
    rather than assigned a made-up function label.
    """
    record = yaml.safe_load(metadata_path.read_text())
    rows = []
    for gene in record.get("genes", []):
        if not gene.get("present", True):
            continue
        role = gene.get("role")
        if not role:
            continue
        rows.append((gene["name"], role))
    return rows


def draw_synteny(
    record_ids: list[str],
    db_root: Path,
    out_path: Path,
    runner: Callable = subprocess.run,
) -> Path:
    """Render a multi-locus synteny/comparison diagram across 2+ curated records via clinker.

    Resolves each `record_id` to its real `locus.gbk` under `db_root`,
    generates clinker's real `--gene_functions`/`--colour_map` CSVs directly
    from every present gene's curated `role` field (see this module's
    docstring for the exact, verified column semantics and the gene-id
    convention clinker itself uses), and shells out to the real `clinker`
    CLI with `-p out_path`.

    Raises `SyntenyRecordError` for any `record_id` that does not resolve to
    exactly one accepted record's `locus.gbk`, and `ValueError` if fewer than
    2 record ids are given (clinker's own purpose is a *comparison*).
    Raises `ClinkerError` if the real `clinker` invocation exits non-zero.
    """
    if len(record_ids) < 2:
        raise ValueError(f"draw_synteny needs 2 or more record ids, got {record_ids!r}")

    db_root = Path(db_root)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    gbk_paths: list[Path] = []
    function_rows: list[tuple[str, str]] = []
    for record_id in record_ids:
        gbk_path, metadata_path = _resolve_record(record_id, db_root)
        gbk_paths.append(gbk_path)
        function_rows.extend(_gene_function_rows(metadata_path))

    roles_seen = sorted({role for _name, role in function_rows})

    gene_functions_path = out_path.parent / f"{out_path.stem}.gene_functions.csv"
    colour_map_path = out_path.parent / f"{out_path.stem}.colour_map.csv"

    with gene_functions_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        for gene_name, role in function_rows:
            writer.writerow([gene_name, role])

    with colour_map_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        for role in roles_seen:
            writer.writerow([role, ROLE_COLORS.get(role, _FALLBACK_COLOR)])

    cmd = [
        "clinker",
        *[str(p) for p in gbk_paths],
        "-gf", str(gene_functions_path),
        "-cm", str(colour_map_path),
        "-p", str(out_path),
        "-f",
    ]
    _run_checked(runner, cmd)
    return out_path
