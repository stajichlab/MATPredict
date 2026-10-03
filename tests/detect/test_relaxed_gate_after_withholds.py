"""The relaxed pass is gated on results that SURVIVE the withhold rules.

Review finding F1 (2026-09-28, results/2026-09-28_fable_review/): the gate
`if not results and relaxed_second_pass` ran BEFORE the modelled-gene bar and
the flank-carried rule. A strict candidate that those rules later withheld
still blocked the relaxed pass. New flank references made 3-gene, 0-modelled
noise clusters reach the strict floor in Rhizomucor pusillus FCH_5_7 and
Absidia glauca; both were withheld, and the real loci -- relaxed-pass calls
before -- were never built.
"""
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _model

ORDER = (
    "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
    "    idiomorph_values: [Plus, Minus]\n    taxonomic_scope: [1]\n    genes:\n"
    + "".join(f"      - {{name: f{i}, role: flanking_variable}}\n" for i in range(1, 7))
    + "      - {name: sexP, role: core_MAT, present_in_idiomorphs: [Plus]}\n"
)
KEY = FamilyKey("P", "MAT")
SPANS = {"sexP": ("c1", 1000, 1600), "f1": ("c1", 3000, 4000)}


def _hit(gene, role, contig, start, end):
    return SearchHit(KEY, gene, role, contig, start, end, "+", 90.0, "rec1", "tblastn_genome")


def _run(tmp_path, *, noise):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(ORDER)
    rec = tmp_path / "P" / "Fam" / "rec1"
    rec.mkdir(parents=True)
    rec.joinpath("metadata.yaml").write_text(
        "record_id: rec1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n")

    def localize(*a, **k):
        hits = [_hit("sexP", "core_MAT", "c1", 1000, 1600),
                _hit("f1", "flanking_variable", "c1", 3000, 4000)]
        if noise:
            # 4 of 7 genes on another contig: clears the strict 0.5 floor,
            # but nothing there can be modelled, so the bar withholds it.
            hits += [_hit(f"f{i}", "flanking_variable", "c2", 1000 * i, 1000 * i + 500)
                     for i in range(2, 6)]
        return hits

    def model(gene_name, method):
        if gene_name not in SPANS:
            return None
        contig, s, e = SPANS[gene_name]
        return _model(gene_name, contig, s, e, family_key=KEY,
                      role="core_MAT" if gene_name == "sexP" else "flanking_variable",
                      method=method)

    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **kw: model(gene_name, "miniprot_refine"),
    )


def test_without_noise_the_relaxed_call_is_reported(tmp_path):
    (r,) = _run(tmp_path, noise=False).results
    assert r.detection_pass == "relaxed"


def test_a_strict_candidate_withheld_by_the_bar_does_not_block_the_relaxed_call(tmp_path):
    outcome = _run(tmp_path, noise=True)
    relaxed = [r for r in outcome.results if r.detection_pass == "relaxed"]
    assert len(relaxed) == 1
    assert relaxed[0].contig == "c1"
    assert relaxed[0].idiomorph == "Plus"
    # the withheld strict noise candidate is still accounted for
    assert outcome.suppressed_unpolished >= 1
