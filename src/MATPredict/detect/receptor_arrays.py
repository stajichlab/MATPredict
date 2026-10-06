"""Pheromone-receptor arrays: report-only grouping of STE3-like loci.

Study: analysis/2026-10-06_agaricomycetes-pr-arrays.md (option 1 and 2). In
1,287 quality-passing Agaricomycete genomes 80% of STE3-like arrays are single
loci and the rest are tandem arrays of 2-10 loci, and the same array was
reported as up to 14 separate calls. This module groups the receptor hits of
every family with a `pheromone_precursor_scan` (the PR family) into arrays and
attaches to each PR call the array it sits in and an `receptor_array_support` flag.

A FLAG, never a gate. Nothing here changes whether a locus is called, its
confidence, its verification label or the call counts; it runs after every
call decision and only adds fields to the report.

Arrays occur for BOTH mating and non-mating receptors. The study found arrays
that hold the mating receptors beside paralogs (Coprinopsis cinerea,
Schizophyllum commune), so array membership does not establish that a locus is
a mating receptor, and `receptor_array_support` says only how much evidence sits behind
the array, not which copy mates.

Definitions (ported from results/2026-10-06_agaricomycetes_pr_arrays/):
- STE3-like locus: the family's non-superseded receptor-gene hits, merged per
  strand where they overlap or lie within `LOCUS_MERGE_GAP_BP` (the study
  merged miniprot alignments per strand; raw reference-protein HSPs of one
  gene can be split by an intron or drawn by several references), kept when
  the hits cover >= 50% of one reference receptor and span <= 8 kb (the
  study's miniprot filter; the raw tblastn hits are mostly weak fragments).
- Array: single linkage over loci on one contig, gap <= `ARRAY_GAP_BP` (50 kb)
  between neighbouring locus ends.
- `supported`: the array has >= 2 loci, OR a pheromone-precursor homology hit
  within the family's CAAX window of a member, OR >= 2 distinct strict-CAAX
  ORFs within that window of a member (4% of arrays by chance in the study,
  against 18% for one). Otherwise `unsupported`.
"""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.caax import CAAX_METHOD, scan_config

#: Single-linkage gap between neighbouring receptor loci (study: 50 kb; 80% of
#: loci are singletons at 50 kb, 78% at 200 kb, 85% at 10 kb).
ARRAY_GAP_BP = 50_000
#: Same-strand receptor hits closer than this are one locus.
LOCUS_MERGE_GAP_BP = 300
#: A locus counts as STE3-like only if its hits cover at least this fraction of
#: one reference receptor, and span at most `MAX_LOCUS_SPAN_BP` (the study kept
#: miniprot alignments covering >= 50% of the query and <= 8 kb). Without this
#: the raw tblastn hits (372 receptor HSPs in Trametes versicolor, 67 merged
#: loci, 6 of them with >= 50% coverage) would make about ten times too many
#: single-locus arrays.
MIN_LOCUS_REF_COVERAGE = 0.5
MAX_LOCUS_SPAN_BP = 8_000

SUPPORTED = "supported"
UNSUPPORTED = "unsupported"

REASON_ARRAY_SIZE = "array_size>=2"
REASON_PRECURSOR = "precursor_homology"
REASON_CAAX_ORFS = "caax_orfs>=2"

RECEPTOR_ARRAYS_NOTE = (
    "Receptor arrays occur for both mating and non-mating STE3-like receptors "
    "(the study found arrays holding paralogs beside the mating receptors in "
    "Coprinopsis cinerea and Schizophyllum commune), so array membership does "
    "not establish that a locus is a mating receptor. receptor_array_support is a flag "
    "on the evidence behind an array; it never changes a call, its confidence "
    "or its verification label."
)

#: Columns a per-locus table (`loci.tsv`) adds for a PR call; None/empty for
#: every other call.
LOCI_ARRAY_COLUMNS = (
    "receptor_array_id", "receptor_array_size", "receptor_array_members", "receptor_array_support", "receptor_array_support_reasons",
)


@dataclass(frozen=True)
class ReceptorArray:
    family: str  # "Phylum:Locus"
    contig: str
    start: int
    end: int
    members: tuple[tuple[int, int, str], ...]  # (start, end, strand), by start
    n_caax_orfs: int
    n_precursor_hits: int
    support: str
    reasons: tuple[str, ...]
    n_calls: int = 0
    window_bp: int = 10_000

    @property
    def receptor_array_id(self) -> str:
        return f"{self.family}:{self.contig}:{self.start}-{self.end}"

    @property
    def size(self) -> int:
        return len(self.members)

    def member_strings(self) -> list[str]:
        return [f"{s}-{e}:{st}" for s, e, st in self.members]

    def as_report(self) -> dict:
        return {
            "receptor_array_id": self.receptor_array_id,
            "family": self.family,
            "contig": self.contig,
            "start": self.start,
            "end": self.end,
            "receptor_array_size": self.size,
            "receptor_array_members": self.member_strings(),
            "strict_caax_orfs": self.n_caax_orfs,
            "precursor_homology_hits": self.n_precursor_hits,
            "receptor_array_support": self.support,
            "receptor_array_support_reasons": list(self.reasons),
            "calls": self.n_calls,
        }


def locus_ref_coverage(hits) -> float:
    """Best fraction of one reference receptor covered by a locus's hits.

    tblastn HSPs are exon-sized, so per reference record the aligned lengths
    of its HSPs are summed (an HSP that mostly overlaps one already counted is
    skipped, so a repeated alignment is not counted twice). Hits with no
    reference length (the proteome fast path: one hit per predicted protein)
    use their own `coverage` percentage; a hit with neither is not filtered.
    """
    best = 0.0
    by_ref: dict[tuple, list] = {}
    for h in hits:
        if h.reference_length_aa and h.align_length_aa:
            by_ref.setdefault((h.reference_record_id, h.reference_length_aa), []).append(h)
        elif h.coverage is not None:
            best = max(best, h.coverage / 100.0)
        else:
            return 1.0
    for (_, ref_len), group in by_ref.items():
        counted: list[tuple[int, int]] = []
        total = 0
        for h in sorted(group, key=lambda x: min(x.start, x.end)):
            lo, hi = sorted((h.start, h.end))
            if any(min(hi, e) - max(lo, s) > 0.5 * (hi - lo) for s, e in counted):
                continue
            counted.append((lo, hi))
            total += h.align_length_aa
        best = max(best, total / ref_len)
    return best


def merge_receptor_loci(
    hits, gap: int = LOCUS_MERGE_GAP_BP, min_coverage: float = MIN_LOCUS_REF_COVERAGE,
    max_span: int = MAX_LOCUS_SPAN_BP,
) -> list[tuple[str, int, int, str]]:
    """(contig, start, end, strand) STE3-like loci from receptor hits.

    Hits are merged per strand where they overlap or lie within `gap`; a
    merged locus is kept when its hits cover at least `min_coverage` of one
    reference receptor and it spans at most `max_span` (pass 0 / a huge span
    to disable).
    """
    by: dict[tuple[str, str], list] = {}
    for h in hits:
        by.setdefault((h.contig, h.strand), []).append(h)
    loci = []

    def close(contig, strand, members, cs, ce):
        if ce - cs + 1 <= max_span and locus_ref_coverage(members) >= min_coverage:
            loci.append((contig, cs, ce, strand))

    for (contig, strand), group in by.items():
        group.sort(key=lambda h: min(h.start, h.end))
        members = [group[0]]
        cs, ce = sorted((group[0].start, group[0].end))
        for h in group[1:]:
            s, e = sorted((h.start, h.end))
            if s - ce <= gap:
                members.append(h)
                ce = max(ce, e)
            else:
                close(contig, strand, members, cs, ce)
                members, cs, ce = [h], s, e
        close(contig, strand, members, cs, ce)
    return sorted(loci)


def group_arrays(loci, gap: int = ARRAY_GAP_BP) -> list[list[tuple[str, int, int, str]]]:
    """Single linkage on one contig: a new array starts where the gap to the
    furthest end so far exceeds `gap` (the study's `make_arrays`)."""
    arrays: list[list] = []
    last_contig, cur_end = None, -1
    for loc in sorted(loci):
        contig, start, end, _ = loc
        if contig != last_contig or start - cur_end > gap:
            arrays.append([loc])
            cur_end = end
        else:
            arrays[-1].append(loc)
            cur_end = max(cur_end, end)
        last_contig = contig
    return arrays


def _near(contig, start, end, hit, window) -> bool:
    lo, hi = sorted((hit.start, hit.end))
    return hit.contig == contig and lo <= end + window and hi >= start - window


def _support(size: int, n_caax: int, n_prec: int) -> tuple[str, tuple[str, ...]]:
    met = []
    if size >= 2:
        met.append(REASON_ARRAY_SIZE)
    if n_prec >= 1:
        met.append(REASON_PRECURSOR)
    if n_caax >= 2:
        met.append(REASON_CAAX_ORFS)
    if met:
        return SUPPORTED, tuple(met)
    return UNSUPPORTED, (
        "single_locus", "no_precursor_homology", f"strict_caax_orfs={n_caax}(<2)",
    )


def build_receptor_arrays(families, hits) -> list[ReceptorArray]:
    """Arrays for every family that has a `pheromone_precursor_scan`.

    `hits` are the run's search hits before clustering, including the
    strict-CAAX scan hits (`caax.CAAX_METHOD`). Read-only.
    """
    out: list[ReceptorArray] = []
    for family in families:
        cfg = scan_config(getattr(family, "pheromone_precursor_scan", None))
        if cfg is None:
            continue
        own = [h for h in hits if h.family_key == family.key]
        receptors = [h for h in own if h.gene_name in cfg["receptor_genes"]
                     and h.method != CAAX_METHOD and h.superseded_by is None]
        precursor_names = {
            g["name"] for g in family.genes
            if g.get("gene_class") == "pheromone_precursor" and g["name"] != cfg["gene"]
        }
        precursors = [h for h in own if h.gene_name in precursor_names
                      and h.method != CAAX_METHOD and h.superseded_by is None]
        caax = [h for h in own if h.method == CAAX_METHOD]
        label = f"{family.key.phylum}:{family.key.locus_name}"
        window = cfg["window_bp"]
        for group in group_arrays(merge_receptor_loci(receptors)):
            contig = group[0][0]
            start = min(loc[1] for loc in group)
            end = max(loc[2] for loc in group)

            def near_member(h):
                return any(_near(contig, s, e, h, window) for _, s, e, _ in group)

            orfs = {(h.contig, min(h.start, h.end), max(h.start, h.end), h.strand)
                    for h in caax if near_member(h)}
            n_prec = sum(1 for h in precursors if near_member(h))
            support, reasons = _support(len(group), len(orfs), n_prec)
            out.append(ReceptorArray(
                family=label, contig=contig, start=start, end=end,
                members=tuple((s, e, st) for _, s, e, st in group),
                n_caax_orfs=len(orfs), n_precursor_hits=n_prec,
                support=support, reasons=reasons, window_bp=window,
            ))
    return out


def _pr_labels(families) -> set[str]:
    return {f"{f.key.phylum}:{f.key.locus_name}" for f in families
            if scan_config(getattr(f, "pheromone_precursor_scan", None)) is not None}


def _is_pr_call(result, labels: set[str]) -> bool:
    if f"{result.family_key.phylum}:{result.family_key.locus_name}" in labels:
        return True
    return any(m.get("family") in labels for m in (result.merged_from or []))


def _best_array(result, arrays: list[ReceptorArray]) -> int | None:
    """Index of the array whose loci overlap the call span most (ties: leftmost)."""
    best, best_n = None, 0
    for i, a in enumerate(arrays):
        if a.contig != result.contig:
            continue
        n = sum(1 for s, e, _ in a.members if s <= result.end and e >= result.start)
        if n > best_n:
            best, best_n = i, n
    return best


def attach_array_support(results, families, arrays: list[ReceptorArray]):
    """(results with a `receptor_array` doc on each PR call, arrays with `n_calls`).

    Only `receptor_array` is set; every other field of every result is
    untouched, so calls, tiers and labels are identical with or without it.
    """
    from dataclasses import replace

    labels = _pr_labels(families)
    arrays = list(arrays)
    calls = [0] * len(arrays)
    new = []
    for r in results:
        if not _is_pr_call(r, labels):
            new.append(r)
            continue
        i = _best_array(r, arrays)
        if i is None:
            doc = {"receptor_array_id": None, "receptor_array_size": None, "receptor_array_members": [],
                   "receptor_array_support": None, "receptor_array_support_reasons": []}
        else:
            calls[i] += 1
            a = arrays[i]
            doc = {"receptor_array_id": a.receptor_array_id, "receptor_array_size": a.size,
                   "receptor_array_members": a.member_strings(), "receptor_array_support": a.support,
                   "receptor_array_support_reasons": list(a.reasons)}
        new.append(replace(r, receptor_array=doc))
    return new, [replace(a, n_calls=n) for a, n in zip(arrays, calls)]


def loci_columns(detected_doc: dict) -> dict:
    """The `LOCI_ARRAY_COLUMNS` cells for one `detected` entry of a report.

    `receptor_array_members` and `receptor_array_support_reasons` are `|`-joined, as
    `genes_found` is in loci.tsv. Empty strings for a non-PR call.
    """
    if "receptor_array_id" not in detected_doc:
        return {c: "" for c in LOCI_ARRAY_COLUMNS}
    return {
        "receptor_array_id": detected_doc["receptor_array_id"] or "",
        "receptor_array_size": "" if detected_doc["receptor_array_size"] is None else detected_doc["receptor_array_size"],
        "receptor_array_members": "|".join(detected_doc.get("receptor_array_members") or []),
        "receptor_array_support": detected_doc["receptor_array_support"] or "",
        "receptor_array_support_reasons": "|".join(detected_doc.get("receptor_array_support_reasons") or []),
    }
