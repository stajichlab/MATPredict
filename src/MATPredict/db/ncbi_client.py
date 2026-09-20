"""NCBI E-utilities client: accession resolution and sequence fetch, cache-backed."""
from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from io import StringIO

from Bio import SeqIO

from MATPredict.db.http_cache import CachedFetcher

_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass(frozen=True)
class CdsStructure:
    """A real, per-exon CDS structure parsed from a GenBank flatfile CDS feature.

    `exons` is in this project's 1-based, fully-closed convention (matching
    `fetch_nucleotide_sequence`'s coordinate contract), listed in transcript order:
    ascending genomic coordinate for a plus-strand feature, descending genomic
    coordinate for a minus-strand feature -- exactly the order `db/validate.py`'s
    `_assemble_transcript` expects and does not itself reorder.
    """

    exons: list[tuple[int, int]]
    strand: str
    codon_start: int
    transl_table: int


def to_schema_exons(cds: CdsStructure) -> list[dict[str, int]]:
    """Convert a `CdsStructure`'s `exons` tuple list into the record schema's shape.

    `CdsStructure.exons` is `list[tuple[int, int]]` (this module's internal, terse
    representation), but `db/_schema/metadata.schema.yaml`'s `genes[].exons` requires
    `list[{start: int, end: int}]` -- a list of dicts. This bridges that gap so a
    curator backfilling a gene's real exon structure from `fetch_cds_structure` can
    assign the result straight onto a record's `gene["exons"]` field without
    hand-rolling the tuple-to-dict conversion each time (previously done ad hoc,
    uncommitted, once per backfill script during this branch's Task 4).
    """
    return [{"start": start, "end": end} for start, end in cds.exons]


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

    def fetch_cds_structure(
        self, accession: str, protein_id: str, region: tuple[int, int] | None = None
    ) -> CdsStructure:
        """Parse a GenBank record's real CDS feature into an exon/codon_start/transl_table
        structure, for populating a curated gene's `exons`/`codon_start`/`transl_table`
        schema fields from the actual deposit rather than a hand-transcribed span -- the
        exact bug class (recording the gene/mRNA span or the outer join() bounds instead
        of the real per-exon structure) that caused two real coordinate errors found and
        fixed earlier in this project's curation work.

        `protein_id` disambiguates which CDS feature to use when an accession carries
        more than one CDS (e.g. a multi-gene MAT locus deposit), matched against each
        CDS feature's own `protein_id` qualifier.

        Coordinate conversion: Biopython's `SeqFeature.location`/`CompoundLocation` uses
        a 0-based start, 1-based-inclusive end (Python-slice-like) convention internally;
        this project uses 1-based, fully-closed. Converting requires adding 1 to each
        part's start only -- verified directly against a real fixture (PV763125.2, a
        genuine multi-exon, minus-strand, 5'-partial, codon_start=2 CDS): Biopython's
        raw `(0, 1172)`-style half-open ints for `877..>1172` came back as
        `(876, 1172)`, and `876 + 1 == 877` matches the real GenBank-text coordinate.

        Exon order: empirically verified against that same real fixture that for a
        minus-strand `complement(join(...))` feature, Biopython's
        `location.parts` already iterates in the REVERSE of the GenBank text listing
        order -- i.e. already in descending-genomic-coordinate (transcript, 5'->3')
        order -- matching this project's exon-list convention (see this module's
        `CdsStructure` docstring and `validate.py`'s `_assemble_transcript`) directly.
        No reversal is performed here: reversing would produce ascending order, which
        is wrong for a minus-strand feature under this project's convention.

        `region`, an optional `(start, end)` pair in this project's 1-based, fully-closed
        convention, requests only that subrange via efetch's `seq_start`/`seq_stop`
        instead of the whole accession. This matters for a whole-genome-assembly-scale
        nuccore accession (an NW_/NC_ RefSeq scaffold or chromosome that is itself a
        `CONTIG`-join "master" record pointing at the real underlying INSDC sequence):
        fetching such an accession's full GenBank flatfile with no range returns only
        its `source`/`CONTIG` lines and zero gene/CDS features, but NCBI's efetch
        renders a real, fully-annotated flatfile on demand for any sub-range of it, with
        coordinates reported relative to that sub-range's own start (i.e. position 1 of
        the response is `region[0]` of the real accession) -- empirically verified
        against NW_006267344.1 (the AbH97_2 HD locus scaffold) and NC_006047.2 (the
        Debaryomyces hansenii CBS767 chromosome E RefSeq record), both "master" records
        for which whole-record efetch returns no CDS at all but a `region`-scoped efetch
        around the gene's own recorded span returns the real CDS with correct exon
        coordinates once the `region[0] - 1` offset below is added back. This method's
        own return value is unaffected: coordinates are converted back into the
        accession's real absolute coordinate space before being returned, so a caller
        never needs to know whether `region` was used.
        """
        params = f"db=nuccore&id={accession}&rettype=gb&retmode=text"
        offset = 0
        if region is not None:
            region_start, region_end = region
            params += f"&seq_start={region_start}&seq_stop={region_end}"
            offset = region_start - 1
        url = self._url("efetch.fcgi", params)
        body = self.fetcher.get(url)
        record = SeqIO.read(StringIO(body), "genbank")
        for feature in record.features:
            if feature.type != "CDS":
                continue
            if feature.qualifiers.get("protein_id", [None])[0] != protein_id:
                continue
            location = feature.location
            exons = [(int(part.start) + 1 + offset, int(part.end) + offset) for part in location.parts]
            strand = "-" if location.strand == -1 else "+"
            codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
            transl_table = int(feature.qualifiers.get("transl_table", ["1"])[0])
            return CdsStructure(
                exons=exons, strand=strand, codon_start=codon_start, transl_table=transl_table
            )
        raise ValueError(f"no CDS with protein_id={protein_id!r} found in {accession}")

    def fetch_taxonomy_lineage(self, taxid: int) -> list[int]:
        """Fetch a taxid's NCBI Taxonomy ancestor lineage as a list of taxids (root-first).

        Uses efetch db=taxonomy, whose XML response carries a `LineageEx` list of
        `{TaxId, ScientificName, Rank}` elements -- one per ancestor, ordered from the
        root of the tree down to (but not including) the queried taxid itself. This is
        the numeric-ancestor-chain data `db.taxonomy.resolve_lineage` (a taxonkit
        subprocess wrapper returning a rank-name string) does not provide, and it is
        what `detect.family_registry.route` needs for lineage-aware scope matching --
        a species/strain taxid's family membership is decided by whether any of its
        ancestor taxids (e.g. a subphylum or class rank) is directly listed in a
        family's `taxonomic_scope`.

        Backed by `self.fetcher`'s on-disk cache (keyed by URL, i.e. by taxid), so
        repeat lookups for the same taxid -- within one process or across separate
        runs -- cost one HTTP round trip total, not one per call.
        """
        url = self._url("efetch.fcgi", f"db=taxonomy&id={taxid}&retmode=xml")
        body = self.fetcher.get(url)
        root = ET.fromstring(body)
        lineage: list[int] = []
        for taxon in root.findall(".//LineageEx/Taxon"):
            tax_id_text = taxon.findtext("TaxId")
            if tax_id_text:
                lineage.append(int(tax_id_text))
        return lineage

    def fetch_taxonomy_phylum(self, taxid: int) -> str | None:
        """Fetch the scientific NAME of `taxid`'s phylum-rank ancestor, or None.

        Reads the SAME `efetch db=taxonomy` document `fetch_taxonomy_lineage`
        reads, at the SAME URL -- `LineageEx` carries `{TaxId, ScientificName,
        Rank}` per ancestor, so the phylum is already present in the response
        that a lineage lookup for this taxid has (or will have) fetched.
        Because `self.fetcher` caches on disk keyed by URL, calling both
        methods for one taxid costs exactly ONE HTTP round trip, not two; that
        is what lets `detect.family_registry.route` add a phylum fallback
        without adding a network call to the routing path.

        Returns None when the document declares no phylum-rank ancestor (an
        unclassified or above-phylum taxid), which callers must treat as "phylum
        unknown" rather than as an error.
        """
        url = self._url("efetch.fcgi", f"db=taxonomy&id={taxid}&retmode=xml")
        root = ET.fromstring(self.fetcher.get(url))
        for taxon in root.findall(".//LineageEx/Taxon"):
            if (taxon.findtext("Rank") or "").strip().lower() == "phylum":
                name = (taxon.findtext("ScientificName") or "").strip()
                return name or None
        return None
