"""NCBI E-utilities client: accession resolution and sequence fetch, cache-backed."""
from __future__ import annotations

import json
from dataclasses import dataclass
from io import StringIO

from Bio import SeqIO

from MATPredict.db.http_cache import CachedFetcher

_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass(frozen=True)
class AccessionStatus:
    """Result of resolving an NCBI accession via esummary."""

    accession: str
    resolved: bool
    resolved_version: str | None
    suppressed: bool


@dataclass
class NcbiClient:
    """Thin wrapper around NCBI E-utilities, backed by an injected CachedFetcher."""

    email: str
    api_key: str | None
    fetcher: CachedFetcher

    def _url(self, path: str, params: str) -> str:
        key_param = f"&api_key={self.api_key}" if self.api_key else ""
        return f"{_EUTILS_BASE}/{path}?{params}&email={self.email}{key_param}"

    def resolve_accession(self, accession: str) -> AccessionStatus:
        """Look up an accession's live/suppressed status and current version via esummary."""
        url = self._url("esummary.fcgi", f"db=nuccore&id={accession}&retmode=json")
        body = self.fetcher.get(url)
        data = json.loads(body)
        uid = data["result"]["uids"][0]
        record = data["result"][uid]
        status = record.get("status", "live")
        return AccessionStatus(
            accession=accession,
            resolved=status == "live",
            resolved_version=record.get("accessionversion") if status == "live" else None,
            suppressed=status == "suppressed",
        )

    def fetch_nucleotide_sequence(
        self, accession: str, start: int, end: int, strand: str | None = None
    ) -> str:
        """Fetch a 1-based, fully-closed nucleotide subrange from a nuccore accession.

        `start`/`end` follow this project's coordinate convention (1-based, fully-closed,
        matching `db/gff_export.write_gff3`'s docstring), which is exactly what NCBI
        efetch's `seq_start`/`seq_stop` params expect -- no offset conversion needed.
        `strand="-"` requests efetch's `strand=2`, which returns the reverse complement
        already oriented 5'->3' along the minus strand, so no local revcomp is needed.
        """
        strand_param = "&strand=2" if strand == "-" else ""
        url = self._url(
            "efetch.fcgi",
            f"db=nuccore&id={accession}&rettype=fasta&retmode=text"
            f"&seq_start={start}&seq_stop={end}{strand_param}",
        )
        body = self.fetcher.get(url)
        record = SeqIO.read(StringIO(body), "fasta-blast")
        return str(record.seq)

    def fetch_protein_sequence(self, accession: str) -> str:
        """Fetch a protein accession's sequence as a plain string (no header, no newlines).

        Uses the "fasta-blast" parser, not "fasta": live NCBI efetch responses can carry
        leading comment lines (e.g. rate-limit/usage notices), which the strict "fasta"
        parser rejects outright. "fasta-blast" tolerates '!'/'#'/';'-prefixed comment lines.
        """
        url = self._url("efetch.fcgi", f"db=protein&id={accession}&rettype=fasta&retmode=text")
        body = self.fetcher.get(url)
        record = SeqIO.read(StringIO(body), "fasta-blast")
        return str(record.seq)
