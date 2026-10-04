"""The genetic code a run translates with.

Curator's ruling, 2026-09-21: "fix CUG-Ser1 code now" -- "learn from taxid
(default is 1 without input)" and "allow also a user to provide a genetic code
override and/or we have an input file format for this".

WHY. `tblastn`, `exonerate` and `miniprot` all accept a translation table
(`-db_gencode`, `--geneticcode`, `-T`) and this pipeline passed NONE of them,
so every genome was translated with table 1. The CUG-Ser1 clade (Serinales,
NCBI genetic code 12, "Alternative Yeast Nuclear") reads CTG as serine, not
leucine, and `MTL` was just scoped to Serinales -- 2,368 genomes in the local
library, 1,154 of them annotated table 12.

Measured cost on the real genes, translating the C. albicans MTL CDS both ways:
MTLa1 2 CTG of 211 aa, MTLa2 2 of 202, MTLalpha1 1 of 194 -- 1-2 residues per
protein, 0.5-1.0%. Small, and fixed anyway because it is systematically wrong
in one direction across a whole clade.

The code is in the SAME cached `efetch db=taxonomy` document the router already
reads for lineage and phylum (`<GCId>`), so deriving it costs no extra network
call.
"""
import xml.etree.ElementTree as ET

import pytest

from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.detect.search import _TBLASTN_OUTFMT  # noqa: F401  (import guard)

DOC = """<?xml version="1.0"?>
<TaxaSet><Taxon>
  <TaxId>5476</TaxId>
  <ScientificName>Candida albicans</ScientificName>
  <GeneticCode><GCId>12</GCId><GCName>Alternative Yeast Nuclear</GCName></GeneticCode>
  <LineageEx>
    <Taxon><TaxId>4890</TaxId><ScientificName>Ascomycota</ScientificName><Rank>phylum</Rank></Taxon>
  </LineageEx>
</Taxon></TaxaSet>"""

NO_CODE = """<?xml version="1.0"?>
<TaxaSet><Taxon><TaxId>1</TaxId><ScientificName>x</ScientificName></Taxon></TaxaSet>"""


class _Fetcher:
    def __init__(self, payload):
        self.payload = payload
        self.calls = []

    def get(self, url):
        self.calls.append(url)
        return self.payload


def _client(payload):
    c = NcbiClient.__new__(NcbiClient)
    c.fetcher = _Fetcher(payload)
    c._url = lambda endpoint, query: f"https://x/{endpoint}?{query}"
    return c


def test_the_genetic_code_is_read_from_the_taxonomy_document():
    assert _client(DOC).fetch_taxonomy_genetic_code(5476) == 12


def test_a_document_with_no_genetic_code_returns_none():
    """None means 'unknown', and the caller falls back to 1. Never guess."""
    assert _client(NO_CODE).fetch_taxonomy_genetic_code(1) is None


def test_it_reuses_the_lineage_document_url():
    """Same URL as the lineage/phylum lookups, so this adds no round trip."""
    c = _client(DOC)
    c.fetch_taxonomy_genetic_code(5476)
    c.fetch_taxonomy_phylum(5476)
    assert len(set(c.fetcher.calls)) == 1


# --- the tools actually receive it -----------------------------------------

def _runner(capture):
    class R:
        returncode = 0
        stdout = ""
        stderr = ""
    def run(cmd, **kw):
        capture.append(cmd)
        return R()
    return run


def test_tblastn_is_given_db_gencode(tmp_path):
    from MATPredict.detect.family_registry import Family, FamilyKey
    from MATPredict.detect.search import search_localize
    key = FamilyKey("P", "MAT")
    fam = Family(key=key, vocabulary_type="enum", idiomorph_values=["a"],
                 idiomorph_pattern=None,
                 genes=[{"name": "g1", "role": "core_MAT"}], taxonomic_scope=[1])
    cmds = []
    search_localize(tmp_path / "g.fa", [fam], tmp_path / "r.faa", {"r": key},
                    runner=_runner(cmds), genetic_code=12)
    tblastn = next(c for c in cmds if c[0] == "tblastn")
    assert "-db_gencode" in tblastn
    assert tblastn[tblastn.index("-db_gencode") + 1] == "12"


def test_tblastn_omits_the_flag_for_the_standard_code(tmp_path):
    """Table 1 is the tools' own default; passing it explicitly would churn
    every existing command line for no behavioural change."""
    from MATPredict.detect.family_registry import Family, FamilyKey
    from MATPredict.detect.search import search_localize
    key = FamilyKey("P", "MAT")
    fam = Family(key=key, vocabulary_type="enum", idiomorph_values=["a"],
                 idiomorph_pattern=None,
                 genes=[{"name": "g1", "role": "core_MAT"}], taxonomic_scope=[1])
    cmds = []
    search_localize(tmp_path / "g.fa", [fam], tmp_path / "r.faa", {"r": key},
                    runner=_runner(cmds), genetic_code=1)
    assert "-db_gencode" not in next(c for c in cmds if c[0] == "tblastn")


def _exonerate_args(tmp_path, code, cmds):
    from MATPredict.detect.family_registry import Family, FamilyKey
    from MATPredict.detect.search import polish_with_exonerate
    key = FamilyKey("P", "MAT")
    fam = Family(key=key, vocabulary_type="enum", idiomorph_values=["a"],
                 idiomorph_pattern=None,
                 genes=[{"name": "g1", "role": "core_MAT"}], taxonomic_scope=[1])
    (tmp_path / "g.fa").write_text(">c1\n" + "ACGT" * 100 + "\n")
    return polish_with_exonerate(tmp_path / "g.fa", fam, "g1", tmp_path / "r.faa", {"r": key},
                                 ("c1", 1, 400), runner=_runner(cmds), genetic_code=code)


TABLE_26 = "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"


def test_exonerate_gets_the_ncbi_string_for_a_code_it_does_not_have(tmp_path):
    """Table 26 (Alaninales, CUG = Ala) is not built into exonerate 2.4.0; by id
    it exits 1, which failed 20 BFD genomes in the v0.6.0 Ascomycota run.
    exonerate takes the table as a 64-letter string instead (TCAG order)."""
    cmds = []
    _exonerate_args(tmp_path, 26, cmds)
    exo = next(c for c in cmds if c[0] == "exonerate")
    assert exo[exo.index("--geneticcode") + 1] == TABLE_26


def test_exonerate_is_skipped_only_for_a_code_with_no_ncbi_table(tmp_path):
    cmds = []
    assert _exonerate_args(tmp_path, 99, cmds) is None
    assert not any(c[0] == "exonerate" for c in cmds)


def test_ncbi_code_strings():
    from MATPredict.detect.search import ncbi_code_string
    assert ncbi_code_string(26) == TABLE_26
    # identical to exonerate's built-in table 1 and table 12 strings
    assert ncbi_code_string(1) == "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    assert ncbi_code_string(12) == "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    assert ncbi_code_string(99) is None


@pytest.mark.parametrize("code", [1, 12, None])
def test_exonerate_still_runs_for_codes_it_has(tmp_path, code):
    cmds = []
    _exonerate_args(tmp_path, code, cmds)
    exo = next(c for c in cmds if c[0] == "exonerate")
    if code in (None, 1):
        assert "--geneticcode" not in exo
    else:
        assert exo[exo.index("--geneticcode") + 1] == str(code)


def test_miniprot_is_given_code_26(tmp_path):
    from MATPredict.detect.family_registry import Family, FamilyKey
    from MATPredict.detect.search import polish_with_miniprot
    key = FamilyKey("P", "MAT")
    fam = Family(key=key, vocabulary_type="enum", idiomorph_values=["a"],
                 idiomorph_pattern=None,
                 genes=[{"name": "g1", "role": "core_MAT"}], taxonomic_scope=[1])
    (tmp_path / "g.fa").write_text(">c1\n" + "ACGT" * 100 + "\n")
    cmds = []
    polish_with_miniprot(tmp_path / "g.fa", fam, "g1", tmp_path / "r.faa", {"r": key},
                         ("c1", 1, 400), runner=_runner(cmds), genetic_code=26)
    mp = next(c for c in cmds if c[0] == "miniprot")
    assert mp[mp.index("-T") + 1] == "26"


def test_the_exonerate_code_table_matches_exonerate_2_4_0():
    from MATPredict.detect.search import EXONERATE_GENETIC_CODES, exonerate_supports_code
    assert EXONERATE_GENETIC_CODES == {1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23}
    assert exonerate_supports_code(None) and exonerate_supports_code(12)
    assert not exonerate_supports_code(26)
