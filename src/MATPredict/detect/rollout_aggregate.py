"""Consolidate a genome-scale detection-rollout batch's per-genome YAML reports
into one summary.

Task 3 (`MATPredict.detect.batch_runner.run_batch`) writes each genome's own
`detection_report.yaml` (via `MATPredict.detect.report.write_detection_report`)
under `out_dir/<taxid>_<accession>/`. This module reads a whole batch's worth
of those YAML files and produces one `RolloutSummary`: a per-family tally of
confidence tiers across every genome, every attempted-but-not-detected family
with its reason, and best-effort taxonomic anomalies.

**Known Task 3 gap this module must not choke on**: `run_batch` creates a
genome's output directory before that genome's `run_pipeline` call runs, so a
genome whose pipeline call fails (decompression failure, `run_pipeline`
exception) leaves an empty `out_dir/<taxid>_<accession>/` directory with no
`detection_report.yaml` inside it. `aggregate_reports` treats a missing or
unparsable report path as "no result for this genome" -- it records a
`GenomeReportError` and continues, it never lets that abort the rest of the
aggregation. This is a real, deliberate design gap from Task 3, not something
this module fixes at the source.

**Anomaly detection is best-effort taxonomy, not a live NCBI call.** The genome
identifier this module works with (`<taxid>_<accession>`, matching the batch
runner's per-genome directory name) already carries a taxid, and this repo
already has an offline lineage mechanism: `MATPredict.db.taxonomy.resolve_lineage`
shells out to `taxonkit reformat` against a local NCBI taxonomy dump, so no
new taxonomy client is introduced here. This is a separate, pre-existing
function from `detect.family_registry.route`'s default resolver
(`MATPredict.db.taxonomy.default_lineage_taxids`, an NCBI-efetch-backed
ancestor-taxid lookup via `NcbiClient`) -- both live in `db/taxonomy.py`, but
they are different mechanisms with different failure modes: a `taxonkit`
outage does not affect family routing, and an NCBI-efetch outage does not
affect this module's anomaly detection.
When `taxonkit` is unavailable or a taxid's lineage can't be resolved, that
genome is simply left out of anomaly detection (its order/class group is
`None`) rather than raising -- anomaly detection is a best-effort signal on
top of the tally, not something the rest of this module depends on.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import yaml

from MATPredict.db.taxonomy import resolve_lineage

# `<taxid>_<accession>` -- the exact directory-naming convention `run_batch`
# uses (see batch_runner.py's `tag = f"{genome.taxid}_{genome.accession}"`).
_GENOME_ID_RE = re.compile(r"^(?P<taxid>\d+)_(?P<accession>.+)$")

_LINEAGE_RANK_RE = re.compile(r"(?P<rank>[a-z])__(?P<name>[^;]*)")


@dataclass(frozen=True)
class NotDetectedEntry:
    """One (genome, family) pair the batch attempted but did not detect."""

    genome: str
    family: str
    reason: str


@dataclass(frozen=True)
class GenomeReportError:
    """A genome directory whose report could not be read (missing file, or a
    file that failed to parse) -- Task 3's known "empty directory on
    pipeline failure" gap, or a genuinely corrupt report."""

    genome: str
    reason: str


@dataclass(frozen=True)
class Anomaly:
    """A genome that reported nothing for `family`, even though `family` was
    detected in at least one other genome (`detected_in`) from the same
    order/class (`taxonomic_group`) as this genome.

    This is a signal worth investigating, not necessarily a bug: real
    biological absence, an assembly gap, or a detection-pipeline miss are all
    possible explanations.
    """

    genome: str
    family: str
    taxonomic_group: str
    detected_in: list[str]


@dataclass
class RolloutSummary:
    """One batch's consolidated detection-rollout result."""

    total_genomes: int
    confidence_tally: dict[str, dict[str, int]] = field(default_factory=dict)
    #: Per family, how many calls of each `locus_class`. Beside the confidence
    #: tally rather than folded into it: they answer different questions --
    #: confidence is how sure the call is, class is WHAT KIND of thing was
    #: found. Added 2026-09-21 with `partial_locus`, which counts as a
    #: detection for recall but must never be merged into a headline locus
    #: count without being visible.
    locus_class_tally: dict[str, dict[str, int]] = field(default_factory=dict)
    not_detected: list[NotDetectedEntry] = field(default_factory=list)
    anomalies: list[Anomaly] = field(default_factory=list)
    genome_errors: list[GenomeReportError] = field(default_factory=list)
    #: Genomes whose report says `routing_mode: not_searched` (curator's ruling
    #: 2026-09-26): no family was searched, so they are neither failures nor
    #: not-detected, and they are left out of anomaly detection.
    not_searched: list[str] = field(default_factory=list)

    def to_doc(self) -> dict:
        """Plain-dict form for YAML serialization (`write_rollout_summary`)."""
        return {
            "total_genomes": self.total_genomes,
            "confidence_tally": self.confidence_tally,
            "locus_class_tally": self.locus_class_tally,
            "not_detected": [
                {"genome": n.genome, "family": n.family, "reason": n.reason}
                for n in self.not_detected
            ],
            "anomalies": [
                {
                    "genome": a.genome,
                    "family": a.family,
                    "taxonomic_group": a.taxonomic_group,
                    "detected_in": a.detected_in,
                }
                for a in self.anomalies
            ],
            "genome_errors": [
                {"genome": g.genome, "reason": g.reason} for g in self.genome_errors
            ],
            "not_searched": list(self.not_searched),
        }


def _genome_id_from_path(report_path: Path) -> str:
    """The genome identifier for `report_path`, taken from its parent
    directory name (`out_dir/<taxid>_<accession>/detection_report.yaml`).
    The report YAML itself carries no taxid/accession field (see
    `report.py`'s `_result_doc`/`write_detection_report`), so the directory
    name is the only place this identifier lives.
    """
    return report_path.parent.name


def _lineage_group(genome_id: str, lineage_resolver: Callable[[int], str | None]) -> str | None:
    """The order/class grouping key for `genome_id`, or `None` if it can't be
    determined (bad genome-id shape, or the resolver couldn't resolve a
    lineage for this taxid). Order is preferred over class (more specific);
    class is the fallback when order is blank -- both are treated as
    "same order/class" per the brief's anomaly definition.
    """
    match = _GENOME_ID_RE.match(genome_id)
    if not match:
        return None
    try:
        taxid = int(match.group("taxid"))
    except ValueError:
        return None
    try:
        lineage = lineage_resolver(taxid)
    except Exception:  # noqa: BLE001 - lineage resolution is best-effort, never fatal
        return None
    if not lineage:
        return None
    ranks = dict(_LINEAGE_RANK_RE.findall(lineage))
    order = ranks.get("o", "").strip()
    klass = ranks.get("c", "").strip()
    if order:
        return f"order:{order}"
    if klass:
        return f"class:{klass}"
    return None


def _default_lineage_resolver(taxid: int) -> str | None:
    """Default lineage resolver: `MATPredict.db.taxonomy.resolve_lineage`
    (taxonkit against a local NCBI taxonomy dump -- no network call). This is
    a separate, pre-existing function in `db/taxonomy.py` from the one
    `detect.family_registry.route` uses by default
    (`default_lineage_taxids`, NCBI-efetch-backed); the two are independent
    mechanisms, not a shared one.
    """
    result = resolve_lineage(taxid)
    if not result.is_current:
        return None
    return result.lineage


def aggregate_reports(
    report_paths: list[Path],
    lineage_resolver: Callable[[int], str | None] = _default_lineage_resolver,
) -> RolloutSummary:
    """Read every YAML report in `report_paths` and consolidate them into one
    `RolloutSummary`.

    `report_paths` is expected to be built by the caller (e.g. the
    `rollout-summary` CLI subcommand) from the batch's `out_dir`, one path per
    genome directory (`out_dir/<taxid>_<accession>/detection_report.yaml`),
    whether or not that path actually exists -- a path that doesn't exist, or
    a report file that fails to parse, is recorded as a `GenomeReportError`
    and skipped, never allowed to raise out of this function and abort the
    rest of the batch's aggregation. `total_genomes` is `len(report_paths)`,
    i.e. every genome the batch attempted, including ones with no report to
    show for it.
    """
    confidence_tally: dict[str, dict[str, int]] = {}
    locus_class_tally: dict[str, dict[str, int]] = {}
    not_detected: list[NotDetectedEntry] = []
    genome_errors: list[GenomeReportError] = []
    not_searched: list[str] = []

    # genome -> set of families it detected; family -> set of genomes that
    # detected it. Built while reading, used afterward for anomaly detection.
    detected_by_genome: dict[str, set[str]] = {}
    genomes_by_family: dict[str, set[str]] = {}
    all_genome_ids: list[str] = []

    for report_path in report_paths:
        genome_id = _genome_id_from_path(report_path)
        all_genome_ids.append(genome_id)
        try:
            doc = yaml.safe_load(report_path.read_text())
        except (FileNotFoundError, OSError, yaml.YAMLError) as exc:
            genome_errors.append(GenomeReportError(genome=genome_id, reason=str(exc)))
            continue
        if not doc:
            genome_errors.append(
                GenomeReportError(genome=genome_id, reason="report file is empty")
            )
            continue
        if doc.get("routing_mode") == "not_searched":
            # Never searched: not a failure, not a negative, and it must not
            # enter `detected_by_genome`, or every relative's call would flag
            # it as an anomaly for "missing" a family it never looked for.
            not_searched.append(genome_id)
            continue

        detected_families: set[str] = set()
        for result in doc.get("detected") or []:
            family = result.get("family")
            confidence = result.get("confidence")
            if family is None or confidence is None:
                continue
            detected_families.add(family)
            genomes_by_family.setdefault(family, set()).add(genome_id)
            confidence_tally.setdefault(family, {}).setdefault(confidence, 0)
            confidence_tally[family][confidence] += 1
            # Guarded, not defaulted: a report written before `locus_class`
            # existed must still count toward confidence rather than being
            # silently filed under an invented class name.
            locus_class = result.get("locus_class")
            if locus_class is not None:
                locus_class_tally.setdefault(family, {}).setdefault(locus_class, 0)
                locus_class_tally[family][locus_class] += 1

        for entry in doc.get("not_detected") or []:
            family = entry.get("family")
            reason = entry.get("reason", "")
            if family is None:
                continue
            not_detected.append(
                NotDetectedEntry(genome=genome_id, family=family, reason=reason)
            )

        detected_by_genome[genome_id] = detected_families

    anomalies = _find_anomalies(
        detected_by_genome=detected_by_genome,
        genomes_by_family=genomes_by_family,
        lineage_resolver=lineage_resolver,
    )

    return RolloutSummary(
        total_genomes=len(report_paths),
        confidence_tally=confidence_tally,
        locus_class_tally=locus_class_tally,
        not_detected=not_detected,
        anomalies=anomalies,
        genome_errors=genome_errors,
        not_searched=not_searched,
    )


def _find_anomalies(
    detected_by_genome: dict[str, set[str]],
    genomes_by_family: dict[str, set[str]],
    lineage_resolver: Callable[[int], str | None],
) -> list[Anomaly]:
    """A genome is anomalous for `family` if some other genome in the same
    order/class group detected `family` and this genome did not.

    Lineage lookups are memoized per genome id since several genomes can
    share a group and `lineage_resolver` may shell out to `taxonkit`.
    """
    group_cache: dict[str, str | None] = {}

    def group_of(genome_id: str) -> str | None:
        if genome_id not in group_cache:
            group_cache[genome_id] = _lineage_group(genome_id, lineage_resolver)
        return group_cache[genome_id]

    anomalies: list[Anomaly] = []
    for family, detected_genomes in genomes_by_family.items():
        # Groups (order/class) that detected this family at least once, and
        # by whom, so a flagged genome's report can say who else found it.
        detected_groups: dict[str, list[str]] = {}
        for genome_id in detected_genomes:
            group = group_of(genome_id)
            if group is None:
                continue
            detected_groups.setdefault(group, []).append(genome_id)

        for genome_id in detected_by_genome:
            if genome_id in detected_genomes:
                continue  # this genome DID detect the family -- not anomalous
            group = group_of(genome_id)
            if group is None or group not in detected_groups:
                continue
            anomalies.append(
                Anomaly(
                    genome=genome_id,
                    family=family,
                    taxonomic_group=group,
                    detected_in=sorted(detected_groups[group]),
                )
            )

    anomalies.sort(key=lambda a: (a.genome, a.family))
    return anomalies


def write_rollout_summary(summary: RolloutSummary, out_path: Path) -> None:
    """Write `summary` as YAML to `out_path`.

    Pilot-scale batches (this rollout's target list is 13 genomes) produce a
    summary of a few KB, well below this project's "compress large text
    output by default" threshold (tens-of-MB+) -- no `.gz`/`.zst` handling is
    added here. A future much-larger rollout should revisit this.
    """
    out_path.write_text(yaml.safe_dump(summary.to_doc(), sort_keys=False))
