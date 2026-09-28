"""A pheromone precursor found by a strict-CAAX ORF scan beside a receptor.

Curator's ruling 2026-09-27 (option ii). B/PR loci are often several receptor
copies plus tiny pheromone precursors that tblastn cannot find, so the
2-distinct-gene admission rule never fires (the Russula nobilis record called
nothing, results/2026-09-27_russulaceae_receptor/). A short ORF ending in the
strict CAAX motif within 10 kb of a receptor hit is counted as the second
gene. Measured basis (results/2026-09-27_pheromone_positional/NOTE.md): strict
C[VI][IV][AVMG] within 10 kb flagged 6/9 known mating receptors, 0/25 other
STE3 copies, 2.5% of random windows; the textbook alphabet hit 42%.
"""
import yaml

from MATPredict.detect.caax import (
    CAAX_METHOD, DEFAULT_MOTIF, STATUS_CAAX_ORF, caax_precursor_hits, scan_caax_orfs,
)
from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _model, _write_record

KEY = FamilyKey("P", "PR")
COMP = str.maketrans("ACGT", "TGCA")


def _revcomp(s):
    return s.translate(COMP)[::-1]


def _orf(tail="TGCGTTATTGCT", n_codons=25):
    """ATG + filler codons + a CAAX tail + TAA. Default tail = C V I A."""
    return "ATG" + "GCA" * n_codons + tail + "TAA"


def _genome(tmp_path, seq, name="c1"):
    path = tmp_path / "genome.fa"
    path.write_text(f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    return path


def test_the_default_motif_is_the_strict_one():
    assert DEFAULT_MOTIF == "C[VI][IV][AVMG]"


def test_a_strict_caax_orf_is_found_with_one_based_coordinates():
    orf = _orf()
    seq = "C" * 100 + orf + "C" * 100
    [o] = scan_caax_orfs(seq)
    assert o.strand == "+"
    assert o.start == 101 and o.end == 100 + len(orf)  # Met .. stop, inclusive
    assert o.motif == "CVIA"
    assert o.length_aa == 1 + 25 + 4


def test_a_reverse_strand_orf_is_found():
    orf = _orf()
    seq = "C" * 100 + _revcomp(orf) + "C" * 100
    [o] = scan_caax_orfs(seq)
    assert o.strand == "-"
    assert (o.start, o.end) == (101, 100 + len(orf))


def test_the_textbook_alphabet_alone_is_not_enough():
    """C-L-L-S passes C[AVLIM][AVLIM]X but not the strict motif."""
    assert scan_caax_orfs("C" * 50 + _orf(tail="TGCCTGCTGTCT") + "C" * 50) == []


def test_length_limits_20_to_130_codons():
    too_short = _orf(n_codons=10)          # Met + 10 + CAAX = 15 codons
    ok_long = _orf(n_codons=120)           # 125 codons
    too_long = _orf(n_codons=140)          # Met 145 codons upstream of the stop
    assert scan_caax_orfs("C" * 30 + too_short + "C" * 30) == []
    assert len(scan_caax_orfs("C" * 30 + ok_long + "C" * 30)) == 1
    assert scan_caax_orfs("C" * 30 + too_long + "C" * 30) == []


def _receptor(start, end, contig="c1"):
    return SearchHit(KEY, "pheromone_receptor", "core_MAT", contig, start, end, "+", 55.0,
                     "rec1", "tblastn_genome", coverage=60.0)


CFG = {"gene": "caax_precursor", "receptor_genes": ["pheromone_receptor"],
       "motif": DEFAULT_MOTIF, "window_bp": 10_000, "min_codons": 20, "max_codons": 130}


class _Fam:
    key = KEY
    pheromone_precursor_scan = CFG
    genes = [{"name": "pheromone_receptor", "role": "core_MAT"},
             {"name": "caax_precursor", "role": "core_MAT", "optional": True}]


def test_only_orfs_inside_the_window_are_hits(tmp_path):
    orf = _orf()
    seq = "C" * 5_000 + "A" * 600 + "C" * 3_000 + orf + "C" * 20_000 + orf + "C" * 1_000
    genome = _genome(tmp_path, seq)
    rec = _receptor(5_001, 5_600)
    hits = caax_precursor_hits(genome, [_Fam()], [rec])
    assert len(hits) == 1
    h = hits[0]
    assert h.gene_name == "caax_precursor" and h.method == CAAX_METHOD
    assert h.family_key == KEY and h.role == "core_MAT"
    assert h.start == 8_601


def test_no_scan_without_a_receptor_hit_or_for_a_family_without_the_setting(tmp_path):
    genome = _genome(tmp_path, "C" * 100 + _orf() + "C" * 100)

    class _Off(_Fam):
        pheromone_precursor_scan = None

    assert caax_precursor_hits(genome, [_Fam()], []) == []
    assert caax_precursor_hits(genome, [_Off()], [_receptor(1, 90)]) == []


ORDER = {
    "phylum": "P",
    "loci": [{
        "locus_name": "PR", "vocabulary_type": "pattern", "idiomorph_pattern": "^B[0-9]+$",
        "taxonomic_scope": [1],
        "pheromone_precursor_scan": {"gene": "caax_precursor",
                                     "receptor_genes": ["pheromone_receptor"]},
        "genes": [
            {"name": "pheromone_receptor", "role": "core_MAT"},
            {"name": "caax_precursor", "role": "core_MAT", "optional": True},
        ],
    }],
}


def _order(tmp_path, scan=True):
    doc = yaml.safe_load(yaml.safe_dump(ORDER))
    if not scan:
        del doc["loci"][0]["pheromone_precursor_scan"]
    (tmp_path / "P").mkdir(exist_ok=True)
    (tmp_path / "P" / "order.yml").write_text(yaml.safe_dump(doc))
    _write_record(tmp_path, locus_name="PR")


def test_the_roster_setting_loads_with_defaults(tmp_path):
    _order(tmp_path)
    [fam] = load_all_families(tmp_path)
    cfg = fam.pheromone_precursor_scan
    assert cfg["motif"] == DEFAULT_MOTIF and cfg["window_bp"] == 10_000
    assert (cfg["min_codons"], cfg["max_codons"]) == (20, 130)


def _run(tmp_path, scan=True):
    _order(tmp_path, scan)
    seq = "C" * 1_000 + "A" * 900 + "C" * 2_000 + _orf() + "C" * 5_000
    genome = _genome(tmp_path, seq)
    polished = []

    def polish(**kw):
        polished.append(kw["gene_name"])
        if kw["gene_name"] == "pheromone_receptor":
            return _model("pheromone_receptor", "c1", 1_001, 1_900, family_key=KEY)
        return None

    result = run_pipeline(
        genome_fasta=genome, proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa",
        search_fast_path=lambda *a, **k: [],
        search_localize=lambda *a, **k: [_receptor(1_001, 1_900)],
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    return result, polished


def test_a_receptor_alone_is_not_a_locus_without_the_scan(tmp_path):
    result, _ = _run(tmp_path, scan=False)
    assert [r for r in result.results if r.withheld_reason is None] == []


def test_the_scan_supplies_the_second_gene_and_the_call_is_made(tmp_path):
    result, polished = _run(tmp_path)
    [call] = [r for r in result.results if r.withheld_reason is None]
    assert "caax_precursor" in call.genes_found
    assert call.polished_genes == 2
    [pre] = [e for e in call.gene_evidence if e.gene_name == "caax_precursor"]
    assert pre.method == CAAX_METHOD and pre.status == STATUS_CAAX_ORF
    assert pre.caax_motif == "CVIA" and pre.orf_length_aa == 30 and pre.orf_count == 1
    # The precursor is never sent to the homology polish tools.
    assert "caax_precursor" not in polished


def test_a_scan_precursor_cannot_by_itself_make_a_call_high(tmp_path):
    result, _ = _run(tmp_path)
    [call] = [r for r in result.results if r.withheld_reason is None]
    assert call.confidence != "high"


def test_the_report_carries_the_orf_fields(tmp_path):
    result, _ = _run(tmp_path)
    out = tmp_path / "report.yaml"
    write_detection_report(result, out)
    doc = yaml.safe_load(out.read_text())
    [pre] = [e for e in doc["detected"][0]["gene_evidence"] if e["gene"] == "caax_precursor"]
    assert pre["method"] == CAAX_METHOD and pre["status"] == STATUS_CAAX_ORF
    assert pre["caax_motif"] == "CVIA" and pre["orf_length_aa"] == 30 and pre["orf_count"] == 1
