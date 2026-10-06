"""Pheromone-receptor arrays: report-only grouping of STE3-like loci.

Study: analysis/2026-10-06_agaricomycetes-pr-arrays.md (option 1 and 2). In
1,287 quality-passing Agaricomycete genomes 80% of STE3-like arrays are single
loci and the rest are tandem arrays of 2-10 loci, and the same array was
reported as up to 14 separate calls. This module groups the receptor hits of
every family with a `pheromone_precursor_scan` (the PR family) into arrays and
attaches to each PR call the array it sits in and an `array_support` flag.

A FLAG, never a gate. Nothing here changes whether a locus is called, its
confidence, its verification label or the call counts; it runs after every
call decision and only adds fields to the report.

Arrays occur for BOTH mating and non-mating receptors. The study found arrays
that hold the mating receptors beside paralogs (Coprinopsis cinerea,
Schizophyllum commune), so array membership does not establish that a locus is
a mating receptor, and `array_support` says only how much evidence sits behind
the array, not which copy mates.

Definitions (ported from results/2026-10-06_agaricomycetes_pr_arrays/):
- STE3-like locus: the family's non-superseded receptor-gene hits, merged per
  strand where they overlap or lie within `LOCUS_MERGE_GAP_BP` (the study
  merged miniprot alignments per strand; raw reference-protein HSPs of one
  gene can be split by an intron or drawn by several references).
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

SUPPORTED = "supported"
UNSUPPORTED = "unsupported"

REASON_ARRAY_SIZE = "array_size>=2"
REASON_PRECURSOR = "precursor_homology"
REASON_CAAX_ORFS = "caax_orfs>=2"

RECEPTOR_ARRAYS_NOTE = (
    "Receptor arrays occur for both mating and non-mating STE3-like receptors "
    "(the study found arrays holding paralogs beside the mating receptors in "
    "Coprinopsis cinerea and Schizophyllum commune), so array membership does "
    "not establish that a locus is a mating receptor. array_support is a flag "
    "on the evidence behind an array; it never changes a call, its confidence "
    "or its verification label."
)

#: Columns a per-locus table (`loci.tsv`) adds for a PR call; None/empty for
#: every other call.
LOCI_ARRAY_COLUMNS = (
    "array_id", "array_size", "array_members", "array_support", "array_support_reasons",
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
    def array_id(self) -> str:
        return f"{self.family}:{self.contig}:{self.start}-{self.end}"

    @property
    def size(self) -> int:
        return len(self.members)

    def member_strings(self) -> list[str]:
        return [f"{s}-{e}:{st}" for s, e, st in self.members]

    def as_report(self) -> dict:
        return {
            "array_id": self.array_id,
            "family": self.family,
            "contig": self.contig,
            "start": self.start,
            "end": self.end,
            "array_size": self.size,
            "array_members": self.member_strings(),
            "strict_caax_orfs": self.n_caax_orfs,
            "precursor_homology_hits": self.n_precursor_hits,
            "array_support": self.support,
            "array_support_reasons": list(self.reasons),
            "calls": self.n_calls,
        }


def merge_receptor_loci(hits, gap: int = LOCUS_MERGE_GAP_BP) -> list[tuple[str, int, int, str]]:
    """(contig, start, end, strand) loci from receptor hits, merged per strand."""
    by: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for h in hits:
        lo, hi = sorted((h.start, h.end))
        by.setdefault((h.contig, h.strand), []).append((lo, hi))
    loci = []
    for (contig, strand), spans in by.items():
        spans.sort()
        cs, ce = spans[0]
        for s, e in spans[1:]:
            if s - ce <= gap:
                ce = max(ce, e)
            else:
                loci.append((contig, cs, ce, strand))
                cs, ce = s, e
        loci.append((contig, cs, ce, strand))
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
            doc = {"array_id": None, "array_size": None, "array_members": [],
                   "array_support": None, "array_support_reasons": []}
        else:
            calls[i] += 1
            a = arrays[i]
            doc = {"array_id": a.array_id, "array_size": a.size,
                   "array_members": a.member_strings(), "array_support": a.support,
                   "array_support_reasons": list(a.reasons)}
        new.append(replace(r, receptor_array=doc))
    return new, [replace(a, n_calls=n) for a, n in zip(arrays, calls)]


def loci_columns(detected_doc: dict) -> dict:
    """The `LOCI_ARRAY_COLUMNS` cells for one `detected` entry of a report.

    `array_members` and `array_support_reasons` are `|`-joined, as
    `genes_found` is in loci.tsv. Empty strings for a non-PR call.
    """
    if "array_id" not in detected_doc:
        return {c: "" for c in LOCI_ARRAY_COLUMNS}
    return {
        "array_id": detected_doc["array_id"] or "",
        "array_size": "" if detected_doc["array_size"] is None else detected_doc["array_size"],
        "array_members": "|".join(detected_doc.get("array_members") or []),
        "array_support": detected_doc["array_support"] or "",
        "array_support_reasons": "|".join(detected_doc.get("array_support_reasons") or []),
    }
