from __future__ import annotations

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient

ESUMMARY_LIVE = """{"result": {"uids": ["1"], "1": {"accessionversion": "GCA_000315115.1", "status": "live"}}}"""
ESUMMARY_SUPPRESSED = """{"result": {"uids": ["1"], "1": {"accessionversion": "GCA_999999999.1", "status": "suppressed"}}}"""
EFETCH_FASTA = ">AAB12345.1 sexP protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n"


def _fake_transport(responses):
    def transport(url: str) -> str:
        for key, body in responses.items():
            if key in url:
                return body
        raise AssertionError(f"unexpected url: {url}")
    return transport


def test_resolve_accession_live(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"esummary": ESUMMARY_LIVE}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    status = client.resolve_accession("GCA_000315115.1")
    assert status.resolved is True
    assert status.suppressed is False
    assert status.resolved_version == "GCA_000315115.1"


def test_resolve_accession_suppressed(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"esummary": ESUMMARY_SUPPRESSED}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    status = client.resolve_accession("GCA_999999999.1")
    assert status.resolved is False
    assert status.suppressed is True


def test_fetch_protein_sequence(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"efetch": EFETCH_FASTA}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    seq = client.fetch_protein_sequence("AAB12345.1")
    assert seq.startswith("MKTAYIAKQRQ")
    assert "\n" not in seq


def test_fetch_nucleotide_sequence_plus_strand(tmp_path):
    fetcher = CachedFetcher(
        cache_dir=tmp_path,
        transport=_fake_transport({"db=nuccore": ">EU009461.1:100-109\nACGTACGTAC\n"}),
    )
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    seq = client.fetch_nucleotide_sequence("EU009461.1", 100, 109, "+")
    assert seq == "ACGTACGTAC"


EFETCH_TAXONOMY_XML = """<?xml version="1.0"?>
<TaxaSet>
<Taxon>
<TaxId>5270</TaxId>
<ScientificName>Ustilago maydis</ScientificName>
<LineageEx>
<Taxon><TaxId>4751</TaxId><ScientificName>Fungi</ScientificName><Rank>kingdom</Rank></Taxon>
<Taxon><TaxId>5204</TaxId><ScientificName>Basidiomycota</ScientificName><Rank>phylum</Rank></Taxon>
<Taxon><TaxId>5157</TaxId><ScientificName>Ustilaginomycotina</ScientificName><Rank>subphylum</Rank></Taxon>
<Taxon><TaxId>5259</TaxId><ScientificName>Ustilaginomycetes</ScientificName><Rank>class</Rank></Taxon>
</LineageEx>
</Taxon>
</TaxaSet>
"""


def test_fetch_taxonomy_lineage_parses_lineage_ex_taxids(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"db=taxonomy": EFETCH_TAXONOMY_XML}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    lineage = client.fetch_taxonomy_lineage(5270)
    assert lineage == [4751, 5204, 5157, 5259]


def test_fetch_nucleotide_sequence_minus_strand_requests_strand2(tmp_path):
    captured_urls = []

    def transport(url: str) -> str:
        captured_urls.append(url)
        return ">EU009461.1:100-109 c\nGTACGTACGT\n"

    fetcher = CachedFetcher(cache_dir=tmp_path, transport=transport)
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    seq = client.fetch_nucleotide_sequence("EU009461.1", 100, 109, "-")
    assert seq == "GTACGTACGT"
    assert "strand=2" in captured_urls[0]


# Real GenBank flatfile for PV763125.2 (Ceratocystis fimbriata MAT1-2-1, partial cds),
# fetched live via efetch during this task's implementation and trimmed to the
# source/CDS feature lines plus the real ORIGIN sequence block (1172 bp, matches the
# feature coordinates exactly). This is the genuine multi-exon, minus-strand, 5'-partial,
# codon_start=2 CDS this project's research identified as the real-world case that
# curators have previously hand-transcribed incorrectly (twice, this session).
GENBANK_FIXTURE_TEXT = """LOCUS       PV763125                1172 bp    DNA     linear   PLN 05-AUG-2026
DEFINITION  Ceratocystis fimbriata isolate Cf089 putative mating type 1-2-1
            protein (MAT1-2-1) gene, partial cds.
ACCESSION   PV763125
VERSION     PV763125.2
KEYWORDS    .
SOURCE      Ceratocystis fimbriata
  ORGANISM  Ceratocystis fimbriata
            Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina;
            Sordariomycetes; Hypocreomycetidae; Microascales;
            Ceratocystidaceae; Ceratocystis.
FEATURES             Location/Qualifiers
     source          1..1172
                     /organism="Ceratocystis fimbriata"
                     /mol_type="genomic DNA"
                     /isolate="Cf089"
                     /db_xref="taxon:5158"
     gene            complement(<1..>1172)
                     /gene="MAT1-2-1"
     CDS             complement(join(156..334,386..684,741..810,877..>1172))
                     /gene="MAT1-2-1"
                     /codon_start=2
                     /product="putative mating type 1-2-1 protein"
                     /protein_id="YGD29557.1"
                     /translation="YPLEMNNTSSFGLHTDLNGIFQPDPQANVNMNYIPSLNHFEMNT
                     IAQSGPTPEMNTVAQFQPGLDMSVVSGSNLGVDSDRVSDINSNANASSAADKIQAVIE
                     ANLVLSLKPKSKNFLLHSTTLTLGTEVHLVRDLQQPHRFLIGDKMLFNTHQKSAVSIP
                     GCEDPLWVEVIPRSLIRPAPQVSKKKVEYRVPRPPNAYILYRKDKHRGVKARNPHMDN
                     NDISIWLGERWRFETSKIRDHYQKTATDYKEMFMLTYPDYQYRPRKANQRKRRAKRAA
                     VSAH"
ORIGIN
        1 aaggcagcaa atccttgtaa atcattcggt agaaaatggt aaatacgaaa ctcattatat
       61 tccatatgca aataaacttc cctagtcaaa ttggtatgaa tgtcatccat cggccctagc
      121 gccgctaata agccaggaac tctgcaagta ggtattcaat gtgccgatac cgcagcccgt
      181 ttggcacggc gctttcgctg gttcgctttc cgggggcgat attgatagtc aggatatgtc
      241 aacatgaaca tttccttgta atctgtggcc gtcttttggt aatggtcccg aatcttcgag
      301 gtttcaaatc tccaccgctc gcctagccat attgctagaa ggttgtcagc atttactcgg
      361 tcgtggcagg ctggggaaag cttacaaata tcattattgt ccatatgagg attcctagcc
      421 ttaacgccac gatgtttgtc tttgcgatac aaaatgtagg cattcggggg gcgaggaact
      481 cgatattcaa ccttcttttt cgaaacctgc ggcgccggcc ggatcaggct tcgaggaatg
      541 acttcgaccc acaacgggtc ttcacatcca ggtatggata cagctgattt ctggtgggta
      601 ttaaagagca tcttatcacc aatgaggaac ctatgcggct gttgaaggtc ccggacaaga
      661 tggacctcag taccgagtgt caggcttcag agtgtgttaa ttagttattc tagtatcatg
      721 ggcaataggt cgagactcac gtcgtgctgt gtagcaggaa attctttgat tttggtttca
      781 atgataatac caagttagct tcaataacgg ctatattaac tgagtaaaac ataaattgcg
      841 gaaacaccat atcaaatatt tttcgcaaaa tcgcaccctg gattttgtca gctgcggatg
      901 aagcattagc gttcgaattg atatcggata ctctatcgga gtccacacca agatttgacc
      961 cagaaacaac actcatatcg agtccgggct gaaattgggc aacagtgttc atctcaggag
     1021 tagggccaga ctgagcgata gtgttcattt cgaagtgatt taatgaaggg atgtaattca
     1081 tattgacatt agcttgaggg tcaggttgga agatgccatt caaatcagtg tgtagaccaa
     1141 atgaagaggt attattcatc tcaagcggat ag
//
"""


def test_fetch_cds_structure_parses_multi_exon_partial_minus_strand_codon_start(tmp_path):
    fetcher = CachedFetcher(
        cache_dir=tmp_path,
        transport=_fake_transport({"db=nuccore&id=PV763125.2&rettype=gb": GENBANK_FIXTURE_TEXT}),
    )
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    result = client.fetch_cds_structure("PV763125.2", protein_id="YGD29557.1")
    # Empirically verified against Biopython's actual parse of this real record (not
    # assumed): for a minus-strand complement(join(...)) CDS, Biopython's
    # SeqFeature.location.parts already iterates in REVERSE of the GenBank text listing
    # order -- i.e. already descending-genomic-coordinate / transcript (5'->3') order,
    # matching this project's exon-list convention (see validate.py's
    # _assemble_transcript and test_validate.py's minus-strand fixtures) directly, with
    # NO further reversal needed in this function.
    assert result.exons == [(877, 1172), (741, 810), (386, 684), (156, 334)]
    assert result.strand == "-"
    assert result.codon_start == 2
    assert result.transl_table == 1


def test_fetch_cds_structure_region_offsets_coordinates_back_to_absolute(tmp_path):
    """A whole-genome-assembly-scale master record (e.g. an NW_/NC_ RefSeq CONTIG-join
    record) returns zero features on a whole-record efetch; `region` requests a
    sub-range instead, whose response reports coordinates relative to that sub-range's
    own start. Task 4's real backfill run hit exactly this against NW_006267344.1 and
    NC_006047.2 -- verify the offset-correction math against this fixture: passing
    region=(1000, 2171) should return the fixture's real 156..1172-space coordinates
    shifted up by 999 (region[0] - 1), recovering the true absolute genomic coordinates,
    while still requesting the windowed efetch (seq_start=1000&seq_stop=2171 in the URL).
    """
    captured_urls: list[str] = []

    def transport(url: str) -> str:
        captured_urls.append(url)
        return GENBANK_FIXTURE_TEXT

    fetcher = CachedFetcher(cache_dir=tmp_path, transport=transport)
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    result = client.fetch_cds_structure("NW_FAKE.1", protein_id="YGD29557.1", region=(1000, 2171))
    assert result.exons == [(1876, 2171), (1740, 1809), (1385, 1683), (1155, 1333)]
    assert result.codon_start == 2
    assert "seq_start=1000" in captured_urls[0]
    assert "seq_stop=2171" in captured_urls[0]


def test_fetch_cds_structure_raises_when_protein_id_not_found(tmp_path):
    fetcher = CachedFetcher(
        cache_dir=tmp_path,
        transport=_fake_transport({"db=nuccore&id=PV763125.2&rettype=gb": GENBANK_FIXTURE_TEXT}),
    )
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    try:
        client.fetch_cds_structure("PV763125.2", protein_id="NOT_A_REAL_PROTEIN.1")
        assert False, "expected ValueError"
    except ValueError as exc:
        assert "NOT_A_REAL_PROTEIN.1" in str(exc)
