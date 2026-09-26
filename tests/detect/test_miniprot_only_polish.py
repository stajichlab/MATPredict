"""A roster gene marked `polish: miniprot` is never sent to exonerate.

Curator's ruling, 2026-09-26. Profiled on C. albicans 3153A with the PAP1/OBP1/
PIK1 flanks added to MTL: 32 of 36 s went to 14 exonerate calls on those long
flank proteins, against 0.5 s for miniprot on the same genes. A 60-genome
ablation with the flanks polished by miniprot alone gave the same genotype in
60/60 genomes at about half the runtime. Scoped per gene in order.yml, so no
other family's flanks change.
"""
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.polish import STATUS_SINGLE

from tests.detect.test_pipeline import _model, _tblastn

ORDER = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n"
    "      - {name: mfa1, role: core_MAT}\n"
    "      - {name: pra1, role: core_MAT}\n"
    "      - {name: PIK1, role: flanking_variable, optional: true, polish: miniprot}\n"
)


def _run(tmp_path, exonerate_calls):
    from tests.detect.test_pipeline import _write_order, _write_record
    _write_order(tmp_path, ORDER)
    _write_record(tmp_path)
    spans = {"mfa1": (100, 200), "pra1": (300, 400), "PIK1": (600, 900)}

    def localize(*a, **k):
        return [_tblastn(n, "c1", *s) for n, s in spans.items()]

    def exonerate(*, gene_name, **kw):
        exonerate_calls.append(gene_name)
        return _model(gene_name, "c1", *spans[gene_name])

    def miniprot(*, gene_name, **kw):
        return _model(gene_name, "c1", *spans[gene_name], method="miniprot_refine")

    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize, polish_with_exonerate=exonerate,
        polish_with_miniprot=miniprot,
    )


def test_a_miniprot_only_gene_skips_exonerate(tmp_path):
    calls = []
    _run(tmp_path, calls)
    assert "PIK1" not in calls
    assert {"mfa1", "pra1"} <= set(calls)


def test_the_gene_is_still_modelled_by_miniprot(tmp_path):
    outcome = _run(tmp_path, [])
    [result] = outcome.results
    evidence = {e.gene_name: e for e in result.gene_evidence}
    assert evidence["PIK1"].status == STATUS_SINGLE
    assert evidence["PIK1"].method == "miniprot_refine"
    assert result.polished_genes == 3
