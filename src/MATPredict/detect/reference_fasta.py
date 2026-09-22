"""Concatenate every accepted record's proteins.faa for search.py to use as a diamond/exonerate database."""
from __future__ import annotations

import re
from pathlib import Path

from MATPredict import logger
from MATPredict.detect.family_registry import (
    FamilyKey,
    load_all_families,
    load_record_families,
)

_HEADER_RE = re.compile(
    r"^>(?P<record_id>[^|]+)\|gene_index=(?P<gene_index>\d+)\|name=(?P<name>[^|]+)\|role=(?P<role>.+)$"
)


def build_reference_fasta(
    db_root: Path, out_path: Path, family_keys: set[FamilyKey] | None = None
) -> Path:
    """Rewrite every db/**/proteins.faa header to `record_id|geneN|name` and concatenate.

    Transforms gff_export.write_proteins_fasta's header format
    `>{record_id}|gene_index={gene_index}|name={name}|role={role}`
    into search.py's `_parse_reference_header` expected format
    `>{record_id}|gene{gene_index}|{name}`.

    `family_keys`, when given, restricts the output to the curated records
    belonging to those families -- i.e. to the families `family_registry.route`
    actually routed this genome to. This file is the tblastn QUERY set, and
    tblastn cost scales with it, so a run narrowed to one phylum must not still
    pay to align every other phylum's proteins: measured against the live
    database, the unrestricted set is 181 proteins while a Mucoromycota-only
    run needs 19. Every protein outside the routed families is not merely
    wasted alignment time, it is a source of spurious cross-phylum candidates
    that the localize/polish stages then grind through one at a time.

    Membership is resolved through `family_registry.load_record_families`,
    which reads each record's own `metadata.yaml` (`phylum` directory +
    `mating_type.locus_name`) -- the same index `search.py` attributes hits
    with, so the query set and the hit attribution cannot disagree. A record
    with no resolvable family is excluded when a restriction is in force: it
    could not be attributed to a routed family anyway.

    Omitting `family_keys` (the default) keeps the pre-existing behaviour
    byte-for-byte, including not reading any `metadata.yaml` at all -- the
    batch/rollout callers, which build one shared reference FASTA for many
    genomes with different taxids, depend on that unrestricted form.
    """
    # Always needed now: the alias map is keyed per family, so a protein's
    # canonical name cannot be resolved without knowing which family its record
    # belongs to -- even on the unrestricted path.
    record_families = load_record_families(db_root)
    _families = load_all_families(db_root)
    aliases_by_family = {f.key: f.gene_aliases for f in _families}
    #: Per family, the canonical gene names the curator has taken OFF the search
    #: list (`exclude_from_search: true` in order.yml). The record keeps the
    #: protein -- this only stops it being used as a query. A gene excluded here
    #: never reaches the reference FASTA, so `searchable_genes_by_family` reads
    #: it back as unsearchable and it leaves the denominator and the core
    #: requirement by the existing route, with no new special case.
    excluded_by_family = {
        f.key: {g["name"] for g in f.genes
                if isinstance(g, dict) and g.get("exclude_from_search")}
        for f in _families
    }
    lines: list[str] = []
    #: (record_id, canonical name, sequence) already emitted. Two curated
    #: proteins that are byte-identical AND resolve to one canonical gene are
    #: one query, not two: MFa1/MFa2/MFa3 are three identical 42 aa entries, and
    #: emitting all three made a single ~95 bp genomic ORF count as three genes.
    #: Identity of SEQUENCE is the test -- MFalpha3, one residue different, is a
    #: distinguishable gene and is emitted separately.
    emitted: set[tuple[str, str, str]] = set()
    for faa in sorted(db_root.glob("*/*/*/proteins.faa")):
        # `db/candidates/` holds not-yet-accepted (needs_review) records, which
        # family_registry.load_record_families and pipeline._short_orf_genes
        # both already exclude as non-authoritative. Including them here would
        # let diamond/exonerate hit a candidate's protein, only for
        # search._attribute to silently drop it later (record_id not found in
        # the accepted-records index) -- wasted search cost with no signal.
        if faa.relative_to(db_root).parts[0] == "candidates":
            continue
        text = faa.read_text()
        for chunk in text.split(">")[1:]:
            header, _, seq = chunk.partition("\n")
            m = _HEADER_RE.match(">" + header)
            if not m:
                logger.warning(f"Skipping malformed header in {faa}: >{header}")
                continue
            # Filtered on the header's own `record_id`, not on the containing
            # directory's name: `load_record_families` is keyed by the
            # `record_id` declared inside `metadata.yaml`, and nothing enforces
            # that a record directory is named after it. Using the header value
            # keeps this filter identical to the key `search._attribute` looks a
            # hit up by, so a protein can never be admitted here and then be
            # unattributable, or vice versa.
            if family_keys is not None and record_families.get(m["record_id"]) not in family_keys:
                continue
            sequence = seq.rstrip("\n")
            family_key = record_families.get(m["record_id"])
            canonical = aliases_by_family.get(family_key, {}).get(m["name"], m["name"])
            if canonical in excluded_by_family.get(family_key, ()):
                continue
            key = (m["record_id"], canonical, sequence.replace("\n", ""))
            if key in emitted:
                continue
            emitted.add(key)
            lines.append(f">{m['record_id']}|gene{m['gene_index']}|{canonical}")
            lines.append(sequence)
    out_path.write_text("\n".join(lines) + "\n")
    return out_path


def searchable_genes_by_family(
    reference_fasta: Path, record_families: dict[str, FamilyKey]
) -> dict[FamilyKey, set[str]]:
    """Per family, the gene names this reference FASTA actually holds a protein for.

    Read back from the WRITTEN file rather than recomputed from `db_root`, so
    the scoring denominator and the query set cannot disagree -- the same
    guarantee `build_reference_fasta` makes for hit attribution by filtering
    on the header's own `record_id`. It matters because phylum routing
    narrows this file: a gene can have a reference somewhere in `db/` and
    still be unsearchable in a given run, and the run must be scored on what
    it could have found, not on what exists elsewhere.

    A record with no resolvable family is skipped; `search._attribute` would
    drop its hits too, so it cannot make any gene findable.

    A missing file yields `{}`, which scoring reads as "no information" and so
    keeps the whole roster in the denominator -- the behaviour that predates
    this function. It warns rather than raising: a caller that stubs out the
    search legitimately never writes this file, but in a real run it is always
    written before the search, so its absence means something else has already
    gone wrong and the run should say so without dying on the scoring step.
    """
    if not reference_fasta.exists():
        logger.warning(
            "No reference FASTA at %s: scoring every expected gene as searchable. "
            "A real run always writes this file before searching.", reference_fasta,
        )
        return {}
    searchable: dict[FamilyKey, set[str]] = {}
    for line in reference_fasta.read_text().splitlines():
        if not line.startswith(">"):
            continue
        record_id, _, gene_name = line[1:].split("|")
        family_key = record_families.get(record_id)
        if family_key is None:
            continue
        searchable.setdefault(family_key, set()).add(gene_name)
    return searchable


def redundant_gene_name_groups(db_root: Path) -> list[dict]:
    """Curated proteins that are byte-identical but resolve to DIFFERENT
    canonical gene names -- i.e. redundancy the roster has not collapsed.

    The general form of the 2026-09-21 ruling. Two identical sequences under
    two gene names cannot be told apart by any homology score: identity,
    coverage, e-value and bitscore are all exactly tied by construction. Left
    in place they inflate `fraction_found`, because one genomic hit is counted
    once per name. Measured before the fix: 927 of 1,384 medium-confidence
    calls across 334 Tremellales genomes were a single ~95 bp ORF reported as
    three MFa or three MFalpha genes.

    Returns one entry per offending group, so a test or a curation check can
    fail with the names in hand. An empty list means every identical-sequence
    group in the database already shares one canonical name.

    Compares within a record, not across records: the same gene curated from
    two strains is legitimately two references for one gene name, and that is
    the ordinary case this must not flag.
    """
    import hashlib
    from collections import defaultdict

    record_families = load_record_families(db_root)
    aliases_by_family = {f.key: f.gene_aliases for f in load_all_families(db_root)}
    by_seq: dict[tuple[str, str], set[str]] = defaultdict(set)
    for faa in sorted(db_root.glob("*/*/*/proteins.faa")):
        if faa.relative_to(db_root).parts[0] == "candidates":
            continue
        for chunk in faa.read_text().split(">")[1:]:
            header, _, seq = chunk.partition("\n")
            m = _HEADER_RE.match(">" + header)
            if not m:
                continue
            family_key = record_families.get(m["record_id"])
            canonical = aliases_by_family.get(family_key, {}).get(m["name"], m["name"])
            digest = hashlib.blake2b(
                seq.strip().replace("\n", "").encode(), digest_size=8
            ).hexdigest()
            by_seq[(m["record_id"], digest)].add(canonical)
    return [
        {"record_id": record_id, "sequence_hash": digest, "canonical_names": sorted(names)}
        for (record_id, digest), names in sorted(by_seq.items())
        if len(names) > 1
    ]
