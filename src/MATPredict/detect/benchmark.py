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

import csv
import gzip
from dataclasses import dataclass
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Seq import Seq

from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.seqmatch import score_match
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.genome_acquisition import DEFAULT_LOCAL_MANIFEST_PATH


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
    * `"exact"` -- the rollout genome unambiguously represents the curated
      record's own source organism, by either of two independent checks:
      (a) **accession match**: the rollout genome's own accession is
      literally one of the curated record's own
      `locus.core.segments[].sequence_source.accession` values (or an
      unversioned prefix of one); or (b) **strain match**: the rollout
      genome and the curated record share a taxid, and the rollout genome's
      own strain (looked up from the local BFD acquisition manifest's
      `STRAIN` column, keyed by its `ASMID`/accession -- see
      `_manifest_strain`) case-insensitively equals the curated record's own
      `organism.strain.name` or one of its `culture_collection_ids`. Either
      path is scored.
    * `"ambiguous"` -- the rollout genome and the curated record share a
      taxid (same species) but NEITHER check above succeeds: the accessions
      differ AND (the manifest strain is unknown/unmounted, or it does not
      match the curated record's strain). **Project judgment call
      (documented in Task 5's report):** this is the real, verified outcome
      for every one of this project's pilot rollout genomes checked against
      its curated Coccidioides/Aspergillus records (5501, 199306, 162425,
      746128) -- the BFD manifest's own `STRAIN` column names a DIFFERENT
      strain (WA_211, 2566, SP-2605-48, niveus) than every existing curated
      record for those species (H538.4, RS, RMSCC1040, Silveira, FGSC A4,
      Af293, A1163), so this is a genuine strain mismatch, not merely an
      accession-namespace mismatch that a same-strain check would have
      resolved. This module EXCLUDES an `"ambiguous"` pairing from the
      numeric sensitivity score rather than assume the specific
      isolate/idiomorph also matches: a same-species genome is not
      guaranteed to carry the same idiomorph (MAT1-1 vs MAT1-2) as the one a
      curator specifically confirmed, and this module's job is to report a
      real, trustworthy number or admit it can't, never a plausible-looking
      guess. Callers that want an "ambiguous" pairing surfaced (e.g. for a
      human to manually confirm strain equivalence when the manifest lookup
      itself is unavailable) get it from this dataclass, never silently
      folded into `FamilyBenchmark.sensitivity`.
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


def _manifest_strain(taxid: int, accession: str, manifest_path: Path) -> str | None:
    """The rollout genome's own strain name, looked up from the local BFD
    acquisition manifest (`genome_acquisition.py`'s `DEFAULT_LOCAL_MANIFEST_PATH`,
    the same CSV `_acquire_local` reads there -- reused here rather than a
    second manifest reader). Matched on `NCBI_TAXONID == str(taxid)` and an
    `ASMID` that is either exactly `accession` or starts with `accession +
    "_"` -- the manifest's real `ASMID` values carry a `<accession>_<assembly
    name>` suffix (e.g. `GCA_004115165.2_Cimm211_ragoo`), while a rollout
    genome's own accession (from its `<taxid>_<accession>` directory name)
    is the bare accession only.

    Returns `None`, never raises, when the manifest isn't mounted (this
    project runs outside UCR HPCC in CI/tests), can't be parsed, or has no
    matching row -- exactly the same "best-effort, not fatal" treatment
    `genome_acquisition._acquire_local` already gives this same file, so a
    missing/unreadable manifest degrades to accession-only matching rather
    than raising out of `match_ground_truth`.
    """
    try:
        with manifest_path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                if row.get("NCBI_TAXONID") != str(taxid):
                    continue
                asmid = row.get("ASMID", "")
                if asmid == accession or asmid.startswith(accession + "_"):
                    strain = row.get("STRAIN")
                    return strain.strip() if strain else None
    except (OSError, csv.Error, KeyError):
        return None
    return None


def _strain_matches(rollout_strain: str, doc: dict) -> bool:
    """Case-insensitive equality between `rollout_strain` and the curated
    record's own `organism.strain.name` or any of its
    `organism.strain.culture_collection_ids` entries."""
    organism_strain = doc.get("organism", {}).get("strain") or {}
    candidates = [organism_strain.get("name")] + list(
        organism_strain.get("culture_collection_ids") or []
    )
    rollout_norm = rollout_strain.strip().lower()
    return any(c and c.strip().lower() == rollout_norm for c in candidates)


def match_ground_truth(
    genome_id: str, db_root: Path, manifest_path: Path = DEFAULT_LOCAL_MANIFEST_PATH,
) -> list[GroundTruthMatch]:
    """Every curated record that shares a taxid with rollout genome
    `genome_id`, classified `"exact"` (accession match OR strain match --
    see `GroundTruthMatch`'s docstring for the exact rule) or `"ambiguous"`
    (same taxid, neither check succeeds -- see `GroundTruthMatch`'s
    docstring for the project's judgment call on why this is excluded from
    the numeric score rather than assumed equivalent).

    A genome with no curated record at all for its taxid returns `[]` --
    there is no ground truth to compare against, ambiguous or otherwise.
    """
    parsed = _genome_id_taxid_accession(genome_id)
    if parsed is None:
        return []
    taxid, accession = parsed
    accession_unversioned = accession.split(".")[0]
    rollout_strain = _manifest_strain(taxid, accession, manifest_path)

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
        if rollout_strain is not None and _strain_matches(rollout_strain, doc):
            matches.append(GroundTruthMatch(
                genome_id=genome_id, record_id=record_id, family_key=key,
                status="exact",
                reason=(
                    f"rollout genome strain {rollout_strain!r} (from the BFD acquisition "
                    f"manifest) matches curated record {record_id}'s own strain, despite a "
                    f"different source accession ({accession} vs "
                    f"{sorted(source_accessions)}) -- same-strain, different assembly"
                ),
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
                f"accession(s) {sorted(source_accessions)}, and rollout strain "
                f"{rollout_strain!r} does not match curated strain {strain!r} "
                "(or the manifest strain lookup was unavailable); "
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


def _splice_transcript(record_seq, exons: list[tuple[int, int]], strand: str | None) -> str:
    """Concatenate `exons` (1-based, fully-closed genomic spans, stored in
    ASCENDING GENOMIC COORDINATE order regardless of strand -- see
    `pipeline.py`'s `GeneEvidence.exons` docstring for why this is NOT
    already transcript order) into one 5'->3' transcript.

    For a plus-strand gene, ascending genomic order already IS transcript
    order. For a minus-strand gene, transcript order is DESCENDING genomic
    order, so the exon list is walked in reverse; each individual exon is
    still reverse-complemented on its own before being appended, exactly the
    per-exon-then-reorder discipline `db/validate.py`'s `_assemble_transcript`
    already documents.
    """
    ordered = list(reversed(exons)) if strand == "-" else list(exons)
    parts = []
    for start, end in ordered:
        span = record_seq[start - 1:end]
        if strand == "-":
            span = span.reverse_complement()
        parts.append(str(span))
    return "".join(parts)


def _extract_translated_gene(
    genome_fasta_path: Path,
    contig: str,
    start: int,
    end: int,
    strand: str | None,
    exons: list[tuple[int, int]] | None = None,
) -> str | None:
    """Re-translate the gene at (contig, start, end, strand) directly from the
    rollout genome's own FASTA. `gene_evidence` entries in a detection report
    carry only coordinates/identity/coverage from the pipeline's own internal
    comparison against its search reference (see `report.py`'s `_result_doc`),
    never the detected gene's own translated protein sequence -- so scoring
    against a curated record's real deposited protein means re-deriving it
    here, exactly once, from the genome FASTA and the reported span.

    When `exons` is a non-empty list (the canonical, polished-model case --
    real MAT-locus genes in these families are routinely multi-exon: e.g.
    COX13 has 5 exons, APN2 has 6, in the Onygenales curated records), the
    real exon structure is spliced via `_splice_transcript` before
    translation -- translating the raw genomic span instead would read
    through introns and produce a systematically wrong-low identity/
    coverage against the curated protein even for a perfectly correct
    detection. When `exons` is `None`/empty (a raw, unpolished hit with no
    exon structure available), this falls back to naive single-span
    translation, which is inherently approximate for a real multi-exon gene
    reported that way -- an accepted limitation of the raw-hit fallback
    path, not of this function.

    A trailing 1-2nt remainder (common for a partial/imprecisely-bounded
    CDS) is trimmed to a multiple of 3 before translation, the same
    discipline `db/validate.py`'s `_independent_translation` already applies
    -- Biopython's `translate()` raises rather than silently truncating one.
    `codon_start` is NOT available on `GeneEvidence` (only the curated
    record's own metadata carries it, and only for the CURATED protein, not
    the rollout's detection), so this assumes `codon_start=1` for the
    rollout side; this is a known, documented simplification, not a defect
    to be silently worked around here.

    Coordinates are this project's 1-based, fully-closed convention (matching
    `gff_export.write_gff3`'s and `NcbiClient.fetch_nucleotide_sequence`'s
    contract). Returns `None` if the contig isn't found in the FASTA.
    """
    with _open_fasta_text(genome_fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id != contig:
                continue
            if exons:
                transcript = _splice_transcript(record.seq, exons, strand)
            else:
                span = record.seq[start - 1:end]
                if strand == "-":
                    span = span.reverse_complement()
                transcript = str(span)
            usable_length = len(transcript) - (len(transcript) % 3)
            return _translate_span(transcript[:usable_length])
    return None


def _gene_evidence_by_name(result_doc: dict) -> dict[str, dict]:
    return {e["gene"]: e for e in result_doc.get("gene_evidence", []) or []}


def score_self_consistency(
    report_paths: list[Path],
    db_root: Path,
    genome_fasta_paths: dict[str, Path],
    ncbi: NcbiClient,
    manifest_path: Path = DEFAULT_LOCAL_MANIFEST_PATH,
) -> tuple[list[FamilyBenchmark], list[GroundTruthMatch]]:
    """Score each rollout genome's detection result against its own curated
    record's ground truth, for genomes that unambiguously match one (see
    `match_ground_truth`).

    For every unambiguous (genome, curated record) pair, for every gene the
    curated record marks `present: true`, the gene counts as "found" iff the
    rollout's detection re-translated protein for that gene (looked up by
    name in the genome's `gene_evidence`, re-translated from the reported
    coordinates -- and real exon structure, when the report carries one --
    against `genome_fasta_paths[genome_id]`) scores `"pass"` or `"warn"`
    (never coordinate/exon comparison) against the curated record's own
    deposited protein (fetched live via `ncbi.fetch_protein_sequence` from
    the gene's `protein_accession`).

    A gene with no matching `gene_evidence` entry at all (the pipeline
    genuinely did not report this gene) counts as a real "not found" --
    that IS the signal self-consistency scoring exists to catch. This is
    kept strictly separate from a gene this function simply COULD NOT
    EVALUATE (no genome FASTA supplied for this genome, the contig wasn't
    found in that FASTA, or the live NCBI protein fetch failed): those genes
    are excluded from the sensitivity denominator entirely rather than
    counted as a fabricated "not found" -- attributing an infrastructure gap
    to the pipeline itself would misrepresent what was actually tested. Each
    not-evaluable gene is recorded (gene name + reason) in the returned
    `FamilyBenchmark.note` so it's visible, never silently dropped.

    Returns `(scored, ground_truth_matches)`:
    * `scored` -- one `FamilyBenchmark` per unambiguously-matched (genome,
      record) pair that had at least one EVALUABLE gene, with a REAL
      `sensitivity` (found / (found + not_found), excluding not-evaluable
      genes from both the numerator and denominator).
      `n_reference_after_holdout` is repurposed here to mean "number of
      self-consistency reference genomes this score is built from" (always 1
      per entry, since each entry is one genome vs one record) -- this
      dataclass is shared with the leave-one-out path by design (see this
      module's docstring); the two fields are self-explanatory for either use
      via `note`.
    * `ground_truth_matches` -- every `GroundTruthMatch` considered,
      `"exact"` and `"ambiguous"` alike, so a caller can see what was
      excluded and why, never silently.

    A genome/record pair with zero EVALUABLE genes (every present gene was
    either genuinely absent from the report -- which alone does not block
    scoring, since 0 found out of N is itself a real result -- or, more to
    the point, one whose genome FASTA/NCBI fetch could not be evaluated at
    all) never emits a `FamilyBenchmark`: nothing to divide by, and this
    function reports a real number or nothing, never a fabricated 0/0.
    """
    curated_docs = {doc["record_id"]: doc for _key, _rid, doc, _path in _load_curated_docs(db_root)}

    scored: list[FamilyBenchmark] = []
    all_matches: list[GroundTruthMatch] = []

    for report_path in report_paths:
        genome_id = report_path.parent.name
        matches = match_ground_truth(genome_id, db_root, manifest_path)
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
            not_found = 0
            not_evaluable: list[str] = []

            for gene in present_genes:
                gene_name = gene["name"]
                evidence = evidence_by_name.get(gene_name)
                if evidence is None:
                    # Genuinely no detection-report evidence for this gene --
                    # a real miss, not an infrastructure gap.
                    not_found += 1
                    continue
                if fasta_path is None:
                    not_evaluable.append(
                        f"{gene_name}: no genome FASTA supplied for {genome_id}"
                    )
                    continue
                exons = (
                    [(e["start"], e["end"]) for e in evidence["exons"]]
                    if evidence.get("exons")
                    else None
                )
                rollout_protein = _extract_translated_gene(
                    fasta_path, evidence["contig"], evidence["start"], evidence["end"],
                    evidence.get("strand"), exons,
                )
                if not rollout_protein:
                    not_evaluable.append(
                        f"{gene_name}: contig {evidence['contig']!r} not found in genome "
                        f"FASTA (or extraction/translation failed)"
                    )
                    continue
                protein_accession = gene.get("protein_accession", "")
                accession = protein_accession.split(":", 1)[-1]
                if not accession:
                    not_evaluable.append(
                        f"{gene_name}: curated record has no protein_accession to compare against"
                    )
                    continue
                try:
                    curated_protein = ncbi.fetch_protein_sequence(accession)
                except Exception as exc:  # noqa: BLE001 - any live NCBI/network
                    # failure here must not abort the rest of this genome's
                    # genes, this genome's other families, or any other
                    # genome in the batch -- one flaky efetch call is an
                    # infrastructure hiccup, not evidence the gene is
                    # missing, so it is excluded (not_evaluable), never
                    # counted as a fabricated miss.
                    not_evaluable.append(
                        f"{gene_name}: NCBI protein fetch failed for {accession!r}: {exc}"
                    )
                    continue
                if score_match(rollout_protein, curated_protein).status in ("pass", "warn"):
                    found += 1
                else:
                    not_found += 1

            evaluable_total = found + not_found
            if evaluable_total == 0:
                continue  # nothing evaluable -- would be a fabricated 0/0

            note = (
                f"self-consistency: rollout genome {genome_id} vs curated record "
                f"{match.record_id} -- {found}/{evaluable_total} evaluable genes "
                f"protein-matched (status={match.status})"
            )
            if not_evaluable:
                note += f"; {len(not_evaluable)} gene(s) not evaluable: " + "; ".join(not_evaluable)

            scored.append(FamilyBenchmark(
                family_key=match.family_key,
                n_reference_after_holdout=1,
                sensitivity=found / evaluable_total,
                note=note,
            ))

    return scored, all_matches
