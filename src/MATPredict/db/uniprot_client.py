"""UniProt REST client: accession resolution and sequence fetch, cache-backed."""
from __future__ import annotations

import json
from dataclasses import dataclass
from io import StringIO

from Bio import SeqIO

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import AccessionStatus

_UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


@dataclass
class UniprotClient:
    """Thin wrapper around the UniProt REST API, backed by an injected CachedFetcher."""

    fetcher: CachedFetcher

    def resolve_accession(self, accession: str) -> AccessionStatus:
        """Confirm a UniProt accession resolves (a 200 JSON entry counts as live)."""
        url = f"{_UNIPROT_BASE}/{accession}.json"
        body = self.fetcher.get(url)
        data = json.loads(body)
        resolved = data.get("primaryAccession") == accession
        return AccessionStatus(accession=accession, resolved=resolved, resolved_version=accession if resolved else None, suppressed=not resolved)

    def fetch_protein_sequence(self, accession: str) -> str:
        """Fetch a UniProt accession's sequence as a plain string.

        Uses "fasta-blast", not "fasta" — see MATPredict.db.ncbi_client.NcbiClient
        for why: some upstream FASTA responses carry leading comment lines that the
        strict "fasta" parser rejects outright.
        """
        url = f"{_UNIPROT_BASE}/{accession}.fasta"
        body = self.fetcher.get(url)
        record = SeqIO.read(StringIO(body), "fasta-blast")
        return str(record.seq)
