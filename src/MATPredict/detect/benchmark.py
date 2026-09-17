"""Species/genus-level leave-one-out Sn/Sp benchmark (sub-project 6)."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from MATPredict.detect.family_registry import FamilyKey


@dataclass(frozen=True)
class FamilyBenchmark:
    family_key: FamilyKey
    n_reference_after_holdout: int
    sensitivity: float | None  # None means "n/a -- insufficient data"
    note: str


def _load_records(db_root: Path) -> list[tuple[FamilyKey, str, Path]]:
    """Return (family_key, species, metadata_path) for every accepted record.

    Accepted records live at db/<Phylum>/<Order-or-Family>/<record_id>/metadata.yaml
    (three path components below db_root). Proposed-but-unaccepted candidates live at
    db/candidates/<Phylum>/<record_id>/metadata.yaml -- also three components below
    db_root, so it matches the same glob shape and must be excluded explicitly by name
    rather than relied on to fall out of a validation.status check (the real candidate
    tree has records at every validation status, including "accepted" ones awaiting
    promotion, so status alone can't distinguish the two trees).
    """
    records = []
    for meta_path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        if meta_path.relative_to(db_root).parts[0] == "candidates":
            continue
        doc = yaml.safe_load(meta_path.read_text())
        key = FamilyKey(meta_path.parents[2].name, doc["mating_type"]["locus_name"])
        species = doc["organism"]["species"]
        records.append((key, species, meta_path))
    return records


def species_genus_holdout_sets(db_root: Path) -> list[tuple[str, list[Path]]]:
    """One (species_or_genus_label, held_out_record_paths) per distinct species/genus
    in the curated DB -- groups same-species Plus/Minus idiomorph pairs into a single
    holdout unit so they can't leak."""
    groups: dict[str, list[Path]] = {}
    for _key, species, path in _load_records(db_root):
        groups.setdefault(species, []).append(path)
    return sorted(groups.items())


def run_benchmark(db_root: Path) -> list[FamilyBenchmark]:
    """For each (phylum, locus_name) family: hold out each species/genus group in turn
    and report whether there is enough curated data to support a leave-one-out
    sensitivity estimate.

    Families with <=2 total curated species report sensitivity=None with an
    explanatory note instead of a misleading score. Real leave-one-out recall
    scoring against `run_pipeline` is a documented follow-up (see the note below);
    this task establishes the holdout-grouping and n/a-reporting contract only.
    """
    records = _load_records(db_root)
    by_family: dict[FamilyKey, list[tuple[str, Path]]] = {}
    for key, species, path in records:
        by_family.setdefault(key, []).append((species, path))

    results: list[FamilyBenchmark] = []
    for key, entries in by_family.items():
        species_groups: dict[str, list[Path]] = {}
        for species, path in entries:
            species_groups.setdefault(species, []).append(path)

        n_groups = len(species_groups)
        if n_groups <= 2:
            results.append(FamilyBenchmark(
                family_key=key, n_reference_after_holdout=max(n_groups - 1, 0),
                sensitivity=None,
                note=f"n/a -- insufficient data ({n_groups} species curated for {key.phylum}:{key.locus_name})",
            ))
            continue

        # Real per-species-group leave-one-out recall scoring is wired here in a
        # follow-up once run_pipeline accepts a pre-built, holdout-filtered
        # reference FASTA; this task establishes the grouping and the
        # n/a-reporting contract the spec requires.
        results.append(FamilyBenchmark(
            family_key=key, n_reference_after_holdout=n_groups - 1,
            sensitivity=None,
            note="holdout grouping ready; recall scoring pending pipeline reference-injection support",
        ))
    return results
