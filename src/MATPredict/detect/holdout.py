"""Which curated records to withhold from a run, so recall can be measured on a
genome the pipeline has not been told the answer for.

**The problem this exists to solve.** A curated record is the QUERY the
pipeline searches with. Running detection on the genome that record came from
asks "did we find the thing we told it to look for", which is a self-consistency
check, not a recall estimate -- `benchmark.py` keeps that distinction carefully
and has reported `sensitivity=None` for every family precisely because nothing
could withhold a record. This module is that missing piece.

**Radius is the dial, and it decides what the number means.** Holding out one
record measures robustness to strain and assembly, which is nearly free to pass
because same-species proteins are ~99% identical. Holding out a whole order
measures the thing that actually matters: whether a clade nobody curated can
still be called. Report the radius beside any recall figure -- a recall number
without one is uninterpretable.

    RECORD   drop this record alone          strain / assembly robustness
    SPECIES  drop every record of its species can a sister species carry it
    GENUS    drop every record of its genus   can the family carry it
    FAMILY   drop every record of its family  }  the real question: an
    ORDER    drop every record of its order   }  uncurated clade

Ranks come from each record's own `taxonomy.lineage`, the same
semicolon-delimited `k__;p__;subphylum__;c__;o__;f__;g__;s__` string the
curation pipeline writes, so no second taxonomy source can drift from it.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import yaml

from MATPredict.detect.family_registry import FamilyKey

logger = logging.getLogger(__name__)


class Radius(str, Enum):
    """How wide a net to cast around the held-out record."""

    RECORD = "record"
    SPECIES = "species"
    GENUS = "genus"
    FAMILY = "family"
    ORDER = "order"


#: Radius -> the lineage prefix that defines it. RECORD has none: it is the
#: record itself, not a taxon.
_RANK_PREFIX = {
    Radius.SPECIES: "s__",
    Radius.GENUS: "g__",
    Radius.FAMILY: "f__",
    Radius.ORDER: "o__",
}


@dataclass(frozen=True)
class RecordTaxon:
    record_id: str
    lineage: dict[str, str]
    path: Path

    def rank(self, prefix: str) -> str | None:
        return self.lineage.get(prefix)


def _parse_lineage(lineage: str | None) -> dict[str, str]:
    """`k__Fungi;p__Ascomycota;...;s__Neurospora_crassa` -> {prefix: value}.

    Tolerant by design: a malformed or absent lineage yields {}, and a record
    with {} can only ever be held out at RECORD radius. Silently widening a
    holdout because a lineage failed to parse would overstate the difficulty of
    the test, which is the direction that flatters the result.
    """
    out: dict[str, str] = {}
    for part in (lineage or "").split(";"):
        part = part.strip()
        if "__" not in part:
            continue
        prefix, _, value = part.partition("__")
        if value:
            out[prefix + "__"] = value
    return out


def load_record_taxa(db_root: Path) -> list[RecordTaxon]:
    """Every ACCEPTED record's id and parsed lineage.

    `db/candidates/` is skipped for the same reason the reference FASTA skips
    it: those records are not authoritative and never reach the query set, so
    they can neither be held out nor leak an answer.
    """
    taxa: list[RecordTaxon] = []
    for meta in sorted(db_root.glob("*/*/*/metadata.yaml")):
        if meta.relative_to(db_root).parts[0] == "candidates":
            continue
        try:
            doc = yaml.safe_load(meta.read_text()) or {}
        except Exception as exc:
            logger.warning("skipping unreadable record %s: %s", meta, exc)
            continue
        rid = doc.get("record_id")
        if not rid:
            continue
        taxa.append(RecordTaxon(rid, _parse_lineage((doc.get("taxonomy") or {}).get("lineage")), meta))
    return taxa


def records_to_withhold(
    db_root: Path, record_id: str, radius: Radius | str
) -> frozenset[str]:
    """Every record id to drop when holding out `record_id` at `radius`.

    Always includes `record_id` itself, even when its lineage is unparseable or
    its rank is missing -- a holdout that failed to withhold its own subject
    would report a self-consistency score as if it were recall, which is the
    exact confusion this module exists to prevent.
    """
    radius = Radius(radius)
    taxa = load_record_taxa(db_root)
    by_id = {t.record_id: t for t in taxa}
    if record_id not in by_id:
        raise KeyError(f"no accepted record {record_id!r} under {db_root}")
    drop = {record_id}
    if radius is not Radius.RECORD:
        prefix = _RANK_PREFIX[radius]
        value = by_id[record_id].rank(prefix)
        if value is None:
            logger.warning(
                "record %s has no %s rank in its lineage; holding out the record "
                "alone rather than silently narrowing the test", record_id, prefix)
        else:
            drop |= {t.record_id for t in taxa if t.rank(prefix) == value}
    return frozenset(drop)


def load_record_idiomorphs(db_root: Path) -> dict[str, set[str]]:
    """record_id -> the idiomorph(s) its `mating_type.idiomorphs` declares.

    A `combined` (homothallic) record declares both, and either is a correct
    call on its locus.
    """
    out: dict[str, set[str]] = {}
    for t in load_record_taxa(db_root):
        try:
            doc = yaml.safe_load(t.path.read_text()) or {}
        except Exception:
            continue
        out[t.record_id] = set((doc.get("mating_type") or {}).get("idiomorphs") or [])
    return out


@dataclass(frozen=True)
class HoldoutScore:
    """How one held-out record fared at one radius.

    `status` is one of:

    * `hit` -- a reported locus overlaps the truth span with a correct idiomorph
    * `hit_undetermined` -- it overlaps, idiomorph undetermined
    * `wrong_idiomorph` -- it overlaps, idiomorph wrong
    * `suppressed` -- nothing reported there, but a locus the modelled-gene
      bar WITHHELD overlaps it: a bar loss, not a search failure
    * `no_reference_family` -- the holdout left no record of the family
    * `no_reference_idiomorph` -- records remain, none for this idiomorph;
      idiomorphs are non-homologous, so the target was unfindable
    * `miss` -- a usable reference remained and nothing was found

    The two `no_reference` kinds are assigned only to a record nothing was
    found for: a locus found anyway (say by its flanks) is scored as found.
    """

    status: str
    locus: dict | None = None


#: Statuses that count as FOUND the locus, for recall. `wrong_idiomorph` is
#: found but wrong; it is reported separately, never folded into recall.
FOUND = frozenset({"hit", "hit_undetermined", "wrong_idiomorph"})
#: Statuses excluded from the recall denominator.
NO_REFERENCE = frozenset({"no_reference_family", "no_reference_idiomorph"})


def _overlapping(loci: list[dict], spans: dict[str, list[int]]) -> list[dict]:
    out = []
    for locus in loci:
        span = spans.get(locus.get("contig"))
        if span and not (locus["end"] < span[0] or locus["start"] > span[1]):
            out.append(locus)
    return out


def score_holdout(
    *,
    record_id: str,
    spans: dict[str, list[int]],
    detected: list[dict],
    suppressed: list[dict],
    withheld: frozenset[str],
    record_families: dict[str, FamilyKey],
    record_idiomorphs: dict[str, set[str]],
) -> HoldoutScore:
    """Score one record's run. `detected`/`suppressed` are report entries
    (`contig`, `start`, `end`, `idiomorph`, ...); `spans` maps a contig to the
    record's true [start, end] on it."""
    expected = record_idiomorphs.get(record_id, set())
    found = _overlapping(detected, spans)
    if found:
        for status, test in (
            ("hit", lambda loc: loc.get("idiomorph") in expected),
            ("hit_undetermined", lambda loc: loc.get("idiomorph") in (None, "undetermined")),
        ):
            match = next((loc for loc in found if test(loc)), None)
            if match is not None:
                return HoldoutScore(status, match)
        return HoldoutScore("wrong_idiomorph", found[0])
    held = _overlapping(suppressed, spans)
    if held:
        return HoldoutScore("suppressed", held[0])
    family = record_families.get(record_id)
    remaining = [r for r, f in record_families.items() if f == family and r not in withheld]
    if not remaining:
        return HoldoutScore("no_reference_family")
    if expected and not any(record_idiomorphs.get(r, set()) & expected for r in remaining):
        return HoldoutScore("no_reference_idiomorph")
    return HoldoutScore("miss")
