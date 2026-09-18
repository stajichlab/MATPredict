from __future__ import annotations

from types import SimpleNamespace
from urllib.parse import parse_qs, urlparse

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.uniprot_client import UniprotClient
from MATPredict.db.validate import _independent_translation, validate_record

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


def _transport_for_ranges(nucleotide_by_range: dict[tuple[str, int, int, str], str]):
    """Build a fake nuccore-efetch transport keyed by explicit (accession, start, stop,
    strand) -> sequence entries, following the same fake-transport-function pattern as
    `_fake_transport` above (a real `NcbiClient` + `CachedFetcher` driven by a fake HTTP
    transport, not a hand-rolled NcbiClient double) but returning fixed per-range strings
    instead of slicing a single flanked chromosome -- needed here because each test
    exercises multiple distinct exon ranges/strands on one fake accession.
    """

    def transport(url: str) -> str:
        if "db=nuccore" not in url:
            raise AssertionError(url)
        params = parse_qs(urlparse(url).query)
        accession = params["id"][0]
        start = int(params["seq_start"][0])
        stop = int(params["seq_stop"][0])
        strand = "-" if params.get("strand") == ["2"] else "+"
        sequence = nucleotide_by_range[(accession, start, stop, strand)]
        return f">{accession}:{start}-{stop}\n{sequence}\n"

    return transport


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


def test_independent_translation_assembles_multi_exon_plus_strand(tmp_path):
    # Two exons, plus strand, codon_start=1. Exon 1: "ATGGCC" (Met-Ala partial),
    # exon 2 continues the frame: "TTTTAA" (Phe-stop). Assembled: ATGGCCTTTTAA
    # -> translates to "MAF" (stop dropped).
    gene = {
        "gene_index": 0, "segment_index": 0, "strand": "+",
        "start": 1, "end": 12,  # outer bounds, unused when exons present
        "exons": [{"start": 1, "end": 6}, {"start": 7, "end": 12}],
        "protein_accession": "ncbi_protein:FAKE1.1",
    }
    record = {"locus": {"core": {"segments": [
        {"sequence_source": {"type": "insdc_nucleotide", "accession": "FAKE_ACC.1"}}
    ]}}}
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_transport_for_ranges({
        ("FAKE_ACC.1", 1, 6, "+"): "ATGGCC",
        ("FAKE_ACC.1", 7, 12, "+"): "TTTTAA",
    }))
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)

    result = _independent_translation(record, gene, ncbi)
    assert result == "MAF"


def test_independent_translation_applies_codon_start_once_not_per_exon(tmp_path):
    # Reproduces the real Ceratocystis PV763125 mechanism: exons on the minus strand,
    # codon_start=2 (skip the first base of the ASSEMBLED sequence, not each exon).
    # Each exon is fetched with strand="-", so `fetch_nucleotide_sequence` already
    # returns it reverse-complemented and 5'->3' oriented (no local revcomp here);
    # exons are listed in transcript order and concatenated as-is. Assembled raw
    # (pre-offset) = "TATGG" + "CCTTTTAA" = "TATGGCCTTTTAA" (13 nt). With
    # codon_start=2, translation starts at index 1: "ATGGCCTTTTAA" -> "MAF" (same
    # result as above, proving the offset is applied once to the whole assembled
    # string, not to each exon's start).
    gene = {
        "gene_index": 0, "segment_index": 0, "strand": "-",
        "start": 1, "end": 13,
        "exons": [{"start": 100, "end": 104}, {"start": 90, "end": 97}],
        "codon_start": 2,
        "protein_accession": "ncbi_protein:FAKE2.1",
    }
    record = {"locus": {"core": {"segments": [
        {"sequence_source": {"type": "insdc_nucleotide", "accession": "FAKE_ACC.1"}}
    ]}}}
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_transport_for_ranges({
        ("FAKE_ACC.1", 100, 104, "-"): "TATGG",
        ("FAKE_ACC.1", 90, 97, "-"): "CCTTTTAA",
    }))
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)

    result = _independent_translation(record, gene, ncbi)
    assert result == "MAF"


def test_independent_translation_falls_back_to_single_span_when_no_exons(tmp_path):
    # Existing behavior (Task 1 did not touch this path): a gene with plain
    # start/end/strand and no `exons` key still works exactly as before.
    gene = {
        "gene_index": 0, "segment_index": 0, "strand": "+",
        "start": 1, "end": 6,
        "protein_accession": "ncbi_protein:FAKE3.1",
    }
    record = {"locus": {"core": {"segments": [
        {"sequence_source": {"type": "insdc_nucleotide", "accession": "FAKE_ACC.1"}}
    ]}}}
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_transport_for_ranges({
        ("FAKE_ACC.1", 1, 6, "+"): "ATGGCC",
    }))
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)

    result = _independent_translation(record, gene, ncbi)
    assert result == "MA"


def test_independent_translation_uses_transl_table_12_for_cug_clade_records(tmp_path):
    # CTG under table 1 is Leu (L); under table 12 (Alternative Yeast Nuclear Code)
    # it is Ser (S). This is the exact real-world failure mode the research found
    # for Candida MTL records.
    gene = {
        "gene_index": 0, "segment_index": 0, "strand": "+",
        "start": 1, "end": 9, "transl_table": 12,
        "protein_accession": "ncbi_protein:FAKE4.1",
    }
    record = {"locus": {"core": {"segments": [
        {"sequence_source": {"type": "insdc_nucleotide", "accession": "FAKE_ACC.1"}}
    ]}}}
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_transport_for_ranges({
        ("FAKE_ACC.1", 1, 9, "+"): "ATGCTGTAA",
    }))
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)

    result = _independent_translation(record, gene, ncbi)
    assert result == "MS"  # not "ML" -- proves table 12 was actually used
