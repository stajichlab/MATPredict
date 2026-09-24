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
