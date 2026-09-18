from __future__ import annotations

from types import SimpleNamespace
from urllib.parse import parse_qs, urlparse

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.uniprot_client import UniprotClient
from MATPredict.db.validate import validate_record

# CDS nucleotide span that translates to exactly "MKTAYIAKQRQISFVKSHFSRQ" + a stop codon,
# embedded in a larger flanking "chromosome" sequence so that fetching the WRONG 1-based
# coordinates genuinely returns different (frameshifted) nucleotides, not just a slice of
# the CDS itself -- used to independently re-derive the protein from the record's own
# coordinates rather than trusting the fetched protein.
_CDS = "ATGAAAACCGCCTATATCGCCAAACAACGTCAAATCAGCTTCGTCAAAAGCCACTTCAGCCGTCAATAA"
_FLANK = "N" * 20
FULL_CHROMOSOME = _FLANK + _CDS + _FLANK
CDS_START = len(_FLANK) + 1  # 1-based
CDS_END = len(_FLANK) + len(_CDS)
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQ"

RECORD = {
    "taxonomy": {"taxid": 4837},
    "locus": {
        "coordinate_provenance": "published_explicit",
        "core": {"segments": [{"segment_index": 0, "sequence_source": {"type": "insdc_nucleotide", "accession": "EU009461.1"}}]},
    },
    "genes": [
        {
            "gene_index": 0,
            "name": "sexP",
            "protein_accession": "ncbi_protein:AAB12345.1",
            "present": True,
            "segment_index": 0,
            "start": CDS_START,
            "end": CDS_END,
            "strand": "+",
        },
    ],
}

ESUMMARY_LIVE = '{"result": {"uids": ["1"], "1": {"accessionversion": "EU009461.1", "status": "live"}}}'
EFETCH_FASTA = f">AAB12345.1\n{PROTEIN}\n"


def _fake_transport(url: str) -> str:
    if "esummary" in url:
        return ESUMMARY_LIVE
    if "db=protein" in url:
        return EFETCH_FASTA
    if "db=nuccore" in url:
        # Mimic real NCBI efetch: seq_start/seq_stop are 1-based, fully-closed.
        params = parse_qs(urlparse(url).query)
        start = int(params["seq_start"][0])
        stop = int(params["seq_stop"][0])
        subseq = FULL_CHROMOSOME[start - 1 : stop]
        return f">EU009461.1:{start}-{stop}\n{subseq}\n"
    raise AssertionError(url)


def _fake_taxonomy_runner(cmd, **kwargs):
    return SimpleNamespace(returncode=0, stdout="4837\tk__Fungi;p__Mucoromycota;...;s__Phycomyces_blakesleeanus\n")


def test_validate_record_all_pass(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport)
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    uniprot = UniprotClient(fetcher=fetcher)

    result = validate_record(
        record={**RECORD, "genes": [{**RECORD["genes"][0], "protein_accession": "ncbi_protein:AAB12345.1"}]},
        ncbi=ncbi,
        uniprot=uniprot,
        taxonomy_runner=_fake_taxonomy_runner,
    )

    assert result["accession_resolved"] is True
    assert result["accession_resolved_version"] == "EU009461.1"
    assert result["sequence_match"]["status"] in {"pass", "warn"}
    assert result["sequence_match"]["per_gene"][0]["percent_identity"] == 100.0
    assert result["taxonomy_current"] is True


def test_validate_record_detects_wrong_coordinates(tmp_path):
    """A regression guard for the tautological `score_match(query=fetched, reference=fetched)`
    bug: that version always fetched one sequence and compared it to itself, so it could
    never fail no matter how wrong the record's curated coordinates were. This record's gene
    coordinates are shifted by +1 nucleotide from the true CDS span (a classic
    0-based/1-based off-by-one), which frameshifts the translation, so the independently
    re-derived translation must NOT match the deposited protein, and the check must fail."""
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport)
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    uniprot = UniprotClient(fetcher=fetcher)

    wrong_gene = {**RECORD["genes"][0], "start": RECORD["genes"][0]["start"] + 1}
    result = validate_record(
        record={**RECORD, "genes": [wrong_gene]},
        ncbi=ncbi,
        uniprot=uniprot,
        taxonomy_runner=_fake_taxonomy_runner,
    )

    assert result["sequence_match"]["status"] == "fail"
    assert result["sequence_match"]["per_gene"][0]["percent_identity"] < 90.0


def test_validate_record_skips_accession_check_for_assembly_type(tmp_path):
    """NcbiClient.resolve_accession only queries NCBI's nuccore database, which cannot
    resolve an assembly accession (GCA_/GCF_) -- that lives in a different NCBI database
    this client doesn't implement lookups for yet. Confirm this is skipped, not crashed.

    The gene here also carries no start/end/strand, so there is no independently
    fetchable/translatable coordinate to check against -- sequence_match must report
    not_applicable rather than a false pass."""
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport)
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    uniprot = UniprotClient(fetcher=fetcher)

    record = {
        "taxonomy": {"taxid": 4837},
        "locus": {
            "coordinate_provenance": "curator_derived",
            "core": {"segments": [{"segment_index": 0, "sequence_source": {"type": "assembly", "accession": "GCA_016772295.1"}}]},
        },
        "genes": [{"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1", "present": True}],
    }

    result = validate_record(record=record, ncbi=ncbi, uniprot=uniprot, taxonomy_runner=_fake_taxonomy_runner)

    assert result["accession_resolved"] is None
    assert result["accession_resolved_version"] is None
    assert result["sequence_match"]["status"] == "not_applicable"
    assert result["sequence_match"]["per_gene"][0]["status"] == "not_applicable"


def test_validate_record_skips_coordinate_checks_when_not_available():
    record = {
        "taxonomy": {"taxid": 4837},
        "locus": {"coordinate_provenance": "not_available"},
        "genes": [],
    }
    result = validate_record(record=record, ncbi=None, uniprot=None, taxonomy_runner=_fake_taxonomy_runner)
    assert result["accession_resolved"] is None
    assert result["taxonomy_current"] is True
    # Regression guard: the early-return path for coordinate_provenance="not_available"
    # used to leave the tautological default sequence_match status of "pass" in place,
    # which is itself a false-pass bug in the same family as the main tautology fix.
    assert result["sequence_match"]["status"] == "not_applicable"
