"""Orchestrates taxonomy/accession/sequence-match checks into a validation result block."""
from __future__ import annotations

import subprocess
from typing import Callable

from Bio.Seq import Seq

from MATPredict.db import taxonomy
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.seqmatch import score_match
from MATPredict.db.uniprot_client import UniprotClient

_STATUS_RANK = {"pass": 0, "warn": 1, "fail": 2}


def _translate_cds(nucleotide_sequence: str, table: int = 1) -> str:
    """Translate a coding sequence to protein, stopping at (and dropping) the first stop
    codon, using the given NCBI genetic code table (default 1, standard code; CUG-clade
    Candida MTL records require table 12, the Alternative Yeast Nuclear Code)."""
    protein = str(Seq(nucleotide_sequence).translate(table=table, to_stop=True))
    return protein


def _assemble_transcript(ncbi: NcbiClient, accession: str, gene: dict) -> str | None:
    """Fetch and concatenate a gene's exon sequences in transcript order, or fall back
    to its single start/end span when no `exons` list is recorded.

    `exons` entries must already be listed in transcript order (5'->3'; for a minus-
    strand gene this is descending genomic coordinate order, matching how GenBank's own
    complement(join(...)) syntax orders its components) -- this function does not
    reorder them. Each exon is fetched with the gene's own strand, so
    `fetch_nucleotide_sequence`'s existing `strand=2` reverse-complement handling
    applies per exon exactly as it already does for a single span; no additional local
    reverse-complementation is performed here.
    """
    strand = gene.get("strand")
    exons = gene.get("exons")
    if exons:
        parts = []
        for exon in exons:
            part = ncbi.fetch_nucleotide_sequence(accession, exon["start"], exon["end"], strand)
            if not part:
                return None
            parts.append(part)
        return "".join(parts)

    start, end = gene.get("start"), gene.get("end")
    if start is None or end is None:
        return None
    return ncbi.fetch_nucleotide_sequence(accession, start, end, strand)


def _independent_translation(record: dict, gene: dict, ncbi: NcbiClient | None) -> str | None:
    """Independently re-derive a gene's protein sequence from the record's own recorded
    genomic coordinates (single span, or a real exon/intron structure when `exons` is
    recorded), by fetching from NCBI and translating -- rather than trusting anything
    already stored in the record about its protein.

    `codon_start` (1/2/3, default 1) is applied EXACTLY ONCE, to the first base of the
    fully assembled, strand-oriented, spliced sequence -- never per-exon. This was
    verified against a real multi-exon, minus-strand, 5'-partial GenBank record; GenBank's
    own /codon_start qualifier is defined relative to "the first base of that feature"
    (the whole joined CDS feature), not per exon.

    Returns None (check not applicable) under the same conditions as before: no NCBI
    client, no fetchable segment sequence, or missing coordinate data on the gene.
    """
    if ncbi is None:
        return None
    segment_index = gene.get("segment_index")
    if segment_index is None:
        return None
    if not gene.get("exons") and (gene.get("start") is None or gene.get("end") is None):
        return None

    segments = record["locus"].get("core", {}).get("segments", [])
    if segment_index >= len(segments):
        return None
    segment = segments[segment_index]
    sequence_source = segment.get("sequence_source", {})
    accession = sequence_source.get("accession")
    if not accession or sequence_source.get("type") != "insdc_nucleotide":
        # Assembly-level accessions (GCA_/GCF_) aren't resolvable by NcbiClient yet --
        # same limitation noted above for accession_resolved. Nothing to independently
        # verify against without fetchable genomic sequence.
        return None

    transcript = _assemble_transcript(ncbi, accession, gene)
    if not transcript:
        return None

    codon_start = gene.get("codon_start", 1)
    offset = codon_start - 1
    frame_corrected = transcript[offset:]
    # A 3'-partial CDS's spliced length is often not a multiple of 3; Biopython's
    # translate() with a trailing 1-2nt remainder needs an explicit trim, since it
    # raises rather than silently truncating.
    usable_length = len(frame_corrected) - (len(frame_corrected) % 3)
    frame_corrected = frame_corrected[:usable_length]

    table = gene.get("transl_table", 1)
    return _translate_cds(frame_corrected, table=table)


def _client_for(accession: str, ncbi: NcbiClient | None, uniprot: UniprotClient | None):
    if accession.startswith("ncbi_protein:"):
        return ncbi, accession.split(":", 1)[1]
    if accession.startswith("uniprotkb:"):
        return uniprot, accession.split(":", 1)[1]
    raise ValueError(f"unrecognized accession namespace: {accession}")


def validate_record(
    record: dict,
    ncbi: NcbiClient | None,
    uniprot: UniprotClient | None,
    taxonomy_runner: Callable = subprocess.run,
) -> dict:
    """Run all applicable validation checks for a candidate record, returning a `validation` dict."""
    result: dict = {
        "accession_resolved": None,
        "accession_resolved_date": None,
        "accession_resolved_version": None,
        "sequence_match": {"status": "not_applicable", "per_gene": [], "notes": ""},
        "taxonomy_current": None,
    }

    taxid = record["taxonomy"]["taxid"]
    tax_result = taxonomy.resolve_lineage(taxid, runner=taxonomy_runner)
    result["taxonomy_current"] = tax_result.is_current

    coordinate_provenance = record["locus"]["coordinate_provenance"]
    if coordinate_provenance == "not_available":
        return result

    segments = record["locus"].get("core", {}).get("segments", [])
    if segments:
        sequence_source = segments[0]["sequence_source"]
        accession = sequence_source.get("accession")
        # NcbiClient.resolve_accession only queries NCBI's nuccore database, which
        # resolves insdc_nucleotide accessions (e.g. "EU009461.1") but not assembly
        # accessions (e.g. "GCA_016772295.1", which live in NCBI's separate assembly
        # database and need a different lookup this client doesn't implement yet).
        # Attempting nuccore esummary on a GCA_/GCF_ accession returns an empty uids
        # list, not an error -- silently skip rather than crash or claim a false result.
        if accession and sequence_source.get("type") == "insdc_nucleotide":
            status = ncbi.resolve_accession(accession)
            result["accession_resolved"] = status.resolved
            result["accession_resolved_version"] = status.resolved_version

    per_gene_results = []
    worst_status = "pass"
    any_checked = False
    for gene in record.get("genes", []):
        if not gene.get("present", True) or not gene.get("protein_accession"):
            continue
        client, bare_accession = _client_for(gene["protein_accession"], ncbi, uniprot)
        fetched_sequence = client.fetch_protein_sequence(bare_accession)

        # Independently re-derive the reference translation from the record's OWN
        # recorded genomic coordinates (segment sequence + gene start/end/strand),
        # rather than comparing the fetched protein to itself (a tautology that always
        # passes and verifies nothing -- see git history). This is the real check: does
        # translating the record's claimed CDS span actually produce the deposited
        # protein?
        independent_translation = _independent_translation(record, gene, ncbi)
        if independent_translation is None:
            per_gene_results.append({
                "gene_index": gene["gene_index"],
                "percent_identity": None,
                "coverage": None,
                "status": "not_applicable",
            })
            continue

        any_checked = True
        match = score_match(query=independent_translation, reference=fetched_sequence)
        per_gene_results.append({
            "gene_index": gene["gene_index"],
            "percent_identity": match.percent_identity,
            "coverage": match.coverage,
            "status": match.status,
        })
        if _STATUS_RANK[match.status] > _STATUS_RANK[worst_status]:
            worst_status = match.status

    overall_status = worst_status if any_checked else "not_applicable"
    result["sequence_match"] = {"status": overall_status, "per_gene": per_gene_results, "notes": ""}
    return result
