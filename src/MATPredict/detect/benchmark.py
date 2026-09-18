"""Species/genus-level leave-one-out Sn/Sp benchmark (sub-project 6), plus
real self-consistency sensitivity scoring against the genome-scale detection
rollout (Task 5 of the 2026-09-18 rollout plan).

**Two different sensitivity computations live in this module, deliberately
kept separate:**

* `run_benchmark` -- the pre-existing species/genus leave-one-out holdout
  grouping across the WHOLE curated DB. It reports `sensitivity=None` for
  every family: real leave-one-out recall scoring needs `run_pipeline` to
  accept a holdout-filtered reference set, which it does not yet support
  (see the note inline below). This function is UNCHANGED by Task 5.
* `score_self_consistency` -- Task 5's new, real, non-stubbed sensitivity
  number. It does not do leave-one-out holdout at all: it takes a genome-scale
  detection rollout's own per-genome YAML reports (Task 4's
  `aggregate_reports` input) and checks, for any rollout genome that is (or
  is a documented-equivalent conspecific/strain match for) an EXISTING
  curated record's own source genome, whether the pipeline re-found -- by
  PROTEIN identity, never by coordinate/exon-structure identity -- the same
  genes the curator already confirmed are really there. This is a
  self-consistency check ("did the pipeline re-derive a fact we already
  know"), not a held-out generalization estimate, so it is intentionally a
  separate function returning the same `FamilyBenchmark` shape rather than a
  branch inside `run_benchmark`: the two numbers answer different questions
  and conflating them under one code path would make it easy to misread a
  self-consistency hit as a held-out recall estimate. A caller (e.g. a future
  rollout-report CLI subcommand) is free to run both and present them
  side by side.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Seq import Seq

from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.seqmatch import score_match
from MATPredict.detect.family_registry import FamilyKey


@dataclass(frozen=True)
class FamilyBenchmark:
    family_key: FamilyKey
    n_reference_after_holdout: int
    sensitivity: float | None  # None means "n/a -- insufficient data"
    note: str


@dataclass(frozen=True)
class GroundTruthMatch:
    """One (rollout genome, curated record) pairing `score_self_consistency`
    considered, together with how confidently the genome represents the
    curated record's own source organism.

    `status` is one of:
    * `"exact"` -- the rollout genome's own accession is literally one of the
      curated record's `locus.core.segments[].sequence_source.accession`
      values (or an unversioned prefix of one). Scored.
    * `"ambiguous"` -- the rollout genome and the curated record share a
      taxid (same species) but the rollout genome's accession does not match
      any of the record's source accessions. **Project judgment call
      (documented in Task 5's report):** this happens for every real pilot
      rollout genome checked against this project's curated Coccidioides/
      Aspergillus records (5501, 199306, 162425, 746128) -- the curated
      records were deposited as locus-level GenBank records (or from a
      different whole-genome assembly) years before the rollout's own
      whole-genome assembly existed, so accessions never coincide even
      though genus/species-level identity is not in doubt. This module
      treats that case as `"ambiguous"` and EXCLUDES it from the numeric
      sensitivity score rather than assume the specific isolate/idiomorph
      also matches: a same-species genome is not guaranteed to carry the
      same idiomorph (MAT1-1 vs MAT1-2) or even the same strain as the one a
      curator specifically confirmed, and this module's job is to report a
      real, trustworthy number or admit it can't, never a plausible-looking
      guess. Callers that want an "ambiguous" pairing surfaced (e.g. for a
      human to manually confirm strain equivalence) get it from this
      dataclass, never silently folded into `FamilyBenchmark.sensitivity`.
    """

    genome_id: str
    record_id: str
    family_key: FamilyKey
    status: str  # "exact" | "ambiguous"
    reason: str


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


# ---------------------------------------------------------------------------
# Task 5: self-consistency sensitivity scoring against the detection rollout.
# ---------------------------------------------------------------------------


def _load_curated_docs(db_root: Path) -> list[tuple[FamilyKey, str, dict, Path]]:
    """Like `_load_records`, but keeps the full parsed metadata doc (taxid,
    genes, source accessions) that self-consistency scoring needs and
    `_load_records` deliberately discards for the leave-one-out path."""
    docs = []
    for meta_path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        if meta_path.relative_to(db_root).parts[0] == "candidates":
            continue
        doc = yaml.safe_load(meta_path.read_text())
        key = FamilyKey(meta_path.parents[2].name, doc["mating_type"]["locus_name"])
        docs.append((key, doc["record_id"], doc, meta_path))
    return docs


def _record_source_accessions(doc: dict) -> set[str]:
    """Every `sequence_source.accession` a curated record's own locus segments
    cite -- the record's ground-truth source genome/accession(s)."""
    accessions = set()
    for segment in doc.get("locus", {}).get("core", {}).get("segments", []) or []:
        accession = segment.get("sequence_source", {}).get("accession")
        if accession:
            accessions.add(accession)
    return accessions


def _genome_id_taxid_accession(genome_id: str) -> tuple[int, str] | None:
    """Parse `<taxid>_<accession>` (the batch runner's directory-naming
    convention, matching `rollout_aggregate._GENOME_ID_RE`) back into its
    parts, or `None` if `genome_id` doesn't have that shape."""
    parts = genome_id.split("_", 1)
    if len(parts) != 2:
        return None
    taxid_part, accession = parts
    try:
        taxid = int(taxid_part)
    except ValueError:
        return None
    return taxid, accession


def match_ground_truth(genome_id: str, db_root: Path) -> list[GroundTruthMatch]:
    """Every curated record that shares a taxid with rollout genome
    `genome_id`, classified `"exact"` (rollout accession literally is one of
    the record's own source accessions) or `"ambiguous"` (same taxid, no
    accession in common -- see `GroundTruthMatch`'s docstring for the
    project's judgment call on why this is excluded from the numeric score
    rather than assumed equivalent).

    A genome with no curated record at all for its taxid returns `[]` --
    there is no ground truth to compare against, ambiguous or otherwise.
    """
    parsed = _genome_id_taxid_accession(genome_id)
    if parsed is None:
        return []
    taxid, accession = parsed
    accession_unversioned = accession.split(".")[0]

    matches: list[GroundTruthMatch] = []
    for key, record_id, doc, _meta_path in _load_curated_docs(db_root):
        if doc.get("taxonomy", {}).get("taxid") != taxid:
            continue
        source_accessions = _record_source_accessions(doc)
        source_unversioned = {a.split(".")[0] for a in source_accessions}
        if accession in source_accessions or accession_unversioned in source_unversioned:
            matches.append(GroundTruthMatch(
                genome_id=genome_id, record_id=record_id, family_key=key,
                status="exact",
                reason=f"rollout genome accession {accession} matches record {record_id}'s "
                       f"own source accession(s) {sorted(source_accessions)}",
            ))
            continue
        species = doc.get("organism", {}).get("species", "?")
        strain = doc.get("organism", {}).get("strain", {}).get("name", "?")
        matches.append(GroundTruthMatch(
            genome_id=genome_id, record_id=record_id, family_key=key,
            status="ambiguous",
            reason=(
                f"same taxid ({taxid}, {species}) but rollout genome accession "
                f"{accession} does not match curated record {record_id}'s own source "
                f"accession(s) {sorted(source_accessions)} (curated strain: {strain}); "
                "treating as ambiguous rather than assuming the same strain/idiomorph "
                "-- excluded from the numeric sensitivity score"
            ),
        ))
    return matches


def _open_fasta_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")  # noqa: SIM115 - caller manages the context


def _translate_span(nucleotide_sequence: str) -> str:
    """Translate a coding sequence to protein, stopping at (and dropping) the
    first stop codon -- the same behaviour as `db.validate._translate_cds`,
    reimplemented here (rather than imported from that private helper) since
    this module needs only the standard genetic-code-table-1 case and
    importing a leading-underscore name from another module across a package
    boundary would tie this module to `db.validate`'s internals for a single
    two-line call."""
    return str(Seq(nucleotide_sequence).translate(table=1, to_stop=True))


def _extract_translated_gene(
    genome_fasta_path: Path, contig: str, start: int, end: int, strand: str | None,
) -> str | None:
    """Re-translate the gene at (contig, start, end, strand) directly from the
    rollout genome's own FASTA. `gene_evidence` entries in a detection report
    carry only coordinates/identity/coverage from the pipeline's own internal
    comparison against its search reference (see `report.py`'s `_result_doc`),
    never the detected gene's own translated protein sequence -- so scoring
    against a curated record's real deposited protein means re-deriving it
    here, exactly once, from the genome FASTA and the reported span.

    Coordinates are this project's 1-based, fully-closed convention (matching
    `gff_export.write_gff3`'s and `NcbiClient.fetch_nucleotide_sequence`'s
    contract). Returns `None` if the contig isn't found in the FASTA.
    """
    with _open_fasta_text(genome_fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id != contig:
                continue
            span = record.seq[start - 1:end]
            if strand == "-":
                span = span.reverse_complement()
            return _translate_span(str(span))
    return None


def _gene_evidence_by_name(result_doc: dict) -> dict[str, dict]:
    return {e["gene"]: e for e in result_doc.get("gene_evidence", []) or []}


def score_self_consistency(
    report_paths: list[Path],
    db_root: Path,
    genome_fasta_paths: dict[str, Path],
    ncbi: NcbiClient,
) -> tuple[list[FamilyBenchmark], list[GroundTruthMatch]]:
    """Score each rollout genome's detection result against its own curated
    record's ground truth, for genomes that unambiguously match one (see
    `match_ground_truth`).

    For every unambiguous (genome, curated record) pair, for every gene the
    curated record marks `present: true`, the gene counts as "found" iff the
    rollout's detection re-translated protein for that gene (looked up by
    name in the genome's `gene_evidence`, re-translated from the reported
    coordinates against `genome_fasta_paths[genome_id]`) scores `"pass"` or
    `"warn"` (never coordinate/exon comparison) against the curated record's
    own deposited protein (fetched live via `ncbi.fetch_protein_sequence` from
    the gene's `protein_accession`). A gene with no matching `gene_evidence`
    entry at all counts as "not found" (score 0), same as one whose protein
    fails the match.

    Returns `(scored, ground_truth_matches)`:
    * `scored` -- one `FamilyBenchmark` per unambiguously-matched (genome,
      record) pair with a REAL `sensitivity` (found / total present genes).
      `n_reference_after_holdout` is repurposed here to mean "number of
      self-consistency reference genomes this score is built from" (always 1
      per entry, since each entry is one genome vs one record) -- this
      dataclass is shared with the leave-one-out path by design (see this
      module's docstring); the two fields are self-explanatory for either use
      via `note`.
    * `ground_truth_matches` -- every `GroundTruthMatch` considered,
      `"exact"` and `"ambiguous"` alike, so a caller can see what was
      excluded and why, never silently.

    A genome/record pair with zero `present: true` genes never happens for a
    real curated record (every accepted record has at least one), but if it
    did, no `FamilyBenchmark` is emitted for it (nothing to divide by) --
    it would otherwise report a fabricated 0/0 "sensitivity".
    """
    curated_docs = {doc["record_id"]: doc for _key, _rid, doc, _path in _load_curated_docs(db_root)}

    scored: list[FamilyBenchmark] = []
    all_matches: list[GroundTruthMatch] = []

    for report_path in report_paths:
        genome_id = report_path.parent.name
        matches = match_ground_truth(genome_id, db_root)
        all_matches.extend(matches)
        exact_matches = [m for m in matches if m.status == "exact"]
        if not exact_matches:
            continue

        try:
            report_doc = yaml.safe_load(report_path.read_text())
        except (FileNotFoundError, OSError, yaml.YAMLError):
            continue
        if not report_doc:
            continue

        results_by_family: dict[str, dict] = {}
        for result in report_doc.get("detected") or []:
            family = result.get("family")
            if family is not None:
                results_by_family[family] = result

        fasta_path = genome_fasta_paths.get(genome_id)

        for match in exact_matches:
            record_doc = curated_docs.get(match.record_id)
            if record_doc is None:
                continue
            family_label = f"{match.family_key.phylum}:{match.family_key.locus_name}"
            result_doc = results_by_family.get(family_label)
            evidence_by_name = _gene_evidence_by_name(result_doc) if result_doc else {}

            present_genes = [g for g in record_doc.get("genes", []) if g.get("present")]
            if not present_genes:
                continue

            found = 0
            for gene in present_genes:
                evidence = evidence_by_name.get(gene["name"])
                if evidence is None or fasta_path is None:
                    continue
                rollout_protein = _extract_translated_gene(
                    fasta_path, evidence["contig"], evidence["start"], evidence["end"],
                    evidence.get("strand"),
                )
                if not rollout_protein:
                    continue
                protein_accession = gene.get("protein_accession", "")
                accession = protein_accession.split(":", 1)[-1]
                if not accession:
                    continue
                curated_protein = ncbi.fetch_protein_sequence(accession)
                if score_match(rollout_protein, curated_protein).status in ("pass", "warn"):
                    found += 1

            scored.append(FamilyBenchmark(
                family_key=match.family_key,
                n_reference_after_holdout=1,
                sensitivity=found / len(present_genes),
                note=(
                    f"self-consistency: rollout genome {genome_id} vs curated record "
                    f"{match.record_id} -- {found}/{len(present_genes)} genes protein-matched "
                    f"(status={match.status})"
                ),
            ))

    return scored, all_matches
