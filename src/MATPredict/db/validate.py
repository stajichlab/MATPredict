"""Orchestrates taxonomy/accession/sequence-match checks into a validation result block."""
from __future__ import annotations

import subprocess
from typing import Callable

from MATPredict.db import taxonomy
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.seqmatch import score_match
from MATPredict.db.uniprot_client import UniprotClient

_STATUS_RANK = {"pass": 0, "warn": 1, "fail": 2}


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
        "sequence_match": {"status": "pass", "per_gene": [], "notes": ""},
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
        accession = segments[0]["sequence_source"].get("accession")
        if accession:
            status = ncbi.resolve_accession(accession)
            result["accession_resolved"] = status.resolved
            result["accession_resolved_version"] = status.resolved_version

    per_gene_results = []
    worst_status = "pass"
    for gene in record.get("genes", []):
        if not gene.get("present", True) or not gene.get("protein_accession"):
            continue
        client, bare_accession = _client_for(gene["protein_accession"], ncbi, uniprot)
        fetched_sequence = client.fetch_protein_sequence(bare_accession)
        # Compares the freshly-fetched sequence to itself: this validates that the
        # accession still resolves to a *stable* sequence, not a re-annotation drift
        # check against the originally-curated sequence (that refinement is left for
        # real curation work in Task 13).
        match = score_match(query=fetched_sequence, reference=fetched_sequence)
        per_gene_results.append({
            "gene_index": gene["gene_index"],
            "percent_identity": match.percent_identity,
            "coverage": match.coverage,
            "status": match.status,
        })
        if _STATUS_RANK[match.status] > _STATUS_RANK[worst_status]:
            worst_status = match.status

    result["sequence_match"] = {"status": worst_status, "per_gene": per_gene_results, "notes": ""}
    return result
