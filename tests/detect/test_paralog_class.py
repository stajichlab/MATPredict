"""R4: a non-MAT HMG paralog class in the idiomorph classifier.

Curator's ruling, J. Stajich, 2026-09-29, on results/2026-09-29_sexM_like_paralog/
NOTE.md. "P1" is a 73-aa HMG gene of the Mucor hiemalis/indicus group that is
present in BOTH mating types (an 85% full-length copy in a Plus-only genome),
sits beside tptA/glrA/algA, and scores 77-84 bits as "sexM" -- enough to pass
the MAT-gene gate on flank support. The replay (p1_replay.tsv) built one HMM
from ONE non-held-out BFD copy (GCA_000697295.1, Mucor indicus B7402): the 23
P1 calls score 145-159 against it and 77-84 against sexM; every other called
core protein scores at least 30 bits LOWER on P1 than on its own idiomorph.

Contract:
* a classifier may carry paralog classes, listed in its manifest
  (`paralog_classes`), each built by the build script from a curated source
  file -- never a hand-edited HMM;
* a verdict records the paralog scores, and names `paralog_class` when the
  best paralog score beats the best MAT score by at least `min_margin`;
* such a call is not a MAT locus: it is withheld with reason `paralog_class`
  (in `suppressed_loci`, with its scores), before the MAT-gene gate;
* without paralog classes nothing changes.
"""
import dataclasses
from pathlib import Path

import pytest
import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.search import SearchHit

from tests.detect.test_idiomorph_classifier import CODON, SEXM, SEXP, _hmm
from tests.detect.test_pipeline import _model

KEY = FamilyKey("P", "MAT")
#: the one BFD P1 copy the replay trained on (GCA_000697295.1, M. indicus B7402)
P1 = "RPLNCFLIYRLEKQQEIVAKCAGANHRDISKIIAKWWKEASEEEKRPFREKARLAKLEHRKMYPGYKYAPRKK"


def _clf(tmp_path, with_paralog=True, min_margin=25.0):
    """An in-memory classifier with real sexP/sexM HMMs and (optionally) P1."""
    import pyhmmer
    from MATPredict.detect.classifier import IdiomorphClassifier
    models, paralogs = {}, {}
    for gene, idiom, seq in (("sexP", "Plus", SEXP), ("sexM", "Minus", SEXM)):
        _hmm(seq, gene, tmp_path / f"{gene}.hmm")
        with pyhmmer.plan7.HMMFile(str(tmp_path / f"{gene}.hmm")) as fh:
            models[idiom] = (gene, fh.read())
    if with_paralog:
        _hmm(P1, "P1", tmp_path / "P1.hmm")
        with pyhmmer.plan7.HMMFile(str(tmp_path / "P1.hmm")) as fh:
            paralogs["P1"] = fh.read()
    return IdiomorphClassifier(directory=tmp_path, models=models, min_margin=min_margin,
                               manifest_sha256="x", paralogs=paralogs)


def test_a_paralog_protein_is_named_by_its_class(tmp_path):
    from MATPredict.detect.classifier import classify
    v = classify(_clf(tmp_path), [P1])
    assert v.paralog_class == "P1"
    assert v.paralog_scores["P1"] - max(v.scores.values()) >= 25.0


def test_a_mat_protein_is_not_a_paralog(tmp_path):
    from MATPredict.detect.classifier import classify
    clf = _clf(tmp_path)
    for seq, idiom in ((SEXP, "Plus"), (SEXM, "Minus")):
        v = classify(clf, [seq])
        assert v.paralog_class is None
        assert v.idiomorph == idiom


def test_without_paralog_classes_the_verdict_is_unchanged(tmp_path):
    from MATPredict.detect.classifier import classify
    v = classify(_clf(tmp_path, with_paralog=False), [P1])
    assert v.paralog_class is None and v.paralog_scores == {}
    assert "paralog_class" not in v.as_report()


def test_the_report_carries_paralog_scores(tmp_path):
    from MATPredict.detect.classifier import classify
    rep = classify(_clf(tmp_path), [P1]).as_report()
    assert rep["paralog_class"] == "P1" and "P1" in rep["paralog_scores"]


def test_combined_verdicts_keep_the_best_paralog_score(tmp_path):
    from MATPredict.detect.classifier import classify, combine_verdicts
    clf = _clf(tmp_path)
    v = combine_verdicts([classify(clf, [P1]), classify(clf, [P1])])
    assert v.paralog_class == "P1"


def test_load_classifier_reads_paralog_classes_from_the_manifest(tmp_path):
    from MATPredict.detect.classifier import ClassifierError, load_classifier
    from MATPredict.detect.family_registry import load_all_families
    from tests.detect.test_idiomorph_classifier import _setup
    _setup(tmp_path)
    cdir = tmp_path / "P" / "classifiers" / "MAT"
    (cdir / "paralogs").mkdir()
    _hmm(P1, "P1", cdir / "paralogs" / "P1.hmm")
    (cdir / "manifest.yaml").write_text(yaml.safe_dump(
        {"family": "P:MAT", "paralog_classes": [{"name": "P1", "hmm": "paralogs/P1.hmm"}]}))
    fam = next(f for f in load_all_families(tmp_path) if f.key == KEY)
    clf = load_classifier(fam.idiomorph_classifier, fam)
    assert set(clf.paralogs) == {"P1"}
    (cdir / "paralogs" / "P1.hmm").unlink()
    from MATPredict.detect import classifier as mod
    mod._CACHE.clear()
    with pytest.raises(ClassifierError):
        load_classifier(fam.idiomorph_classifier, fam)


# --- the gate step -----------------------------------------------------------

def _result(clf_doc, split_locus=None):
    from MATPredict.detect.pipeline import DetectionResult, GeneEvidence
    ev = [GeneEvidence("sexM", "core_MAT", "c1", 1, 100, "+", 38.0, None, "rec1",
                       "exonerate_refine", status="polished_single"),
          GeneEvidence("tptA", "flanking_conserved", "c1", 200, 300, "+", 80.0, None, "rec1",
                       "exonerate_refine", status="polished_single"),
          GeneEvidence("glrA", "flanking_variable", "c1", 400, 500, "+", 80.0, None, "rec1",
                       "exonerate_refine", status="polished_single")]
    return DetectionResult(
        family_key=FamilyKey("Mucoromycota", "MAT"), contig="c1", start=1, end=1000,
        confidence="high", idiomorph="Minus", ambiguous_with=[],
        genes_found=["glrA", "sexM", "tptA"], genes_missing=[], fragmented=False,
        gene_evidence=ev, polished_genes=3, idiomorph_classifier=clf_doc,
        split_locus=split_locus)


P1_DOC = {"method": "hmm", "idiomorph": "Minus", "scores": {"Minus": 78.9, "Plus": 49.2},
          "margin": 29.7, "classifier_input": "model", "min_margin": 25,
          "paralog_scores": {"P1": 158.6}, "paralog_class": "P1"}
SPECS = {FamilyKey("Mucoromycota", "MAT"): {"type": "hmm", "dir": "/x", "min_margin": 25}}


def test_a_paralog_class_call_is_withheld_even_with_flank_support():
    from MATPredict.detect.mat_gene_gate import WITHHELD_PARALOG_CLASS, apply_mat_gene_gate
    kept, withheld = apply_mat_gene_gate([_result(P1_DOC)], SPECS)
    assert kept == []
    (w,) = withheld
    assert w.withheld_reason == WITHHELD_PARALOG_CLASS
    assert w.withheld_detail["paralog_class"] == "P1"
    assert w.withheld_detail["paralog_score"] == 158.6


def test_a_paralog_class_split_locus_call_is_withheld_too():
    from MATPredict.detect.mat_gene_gate import WITHHELD_PARALOG_CLASS, apply_mat_gene_gate
    kept, (w,) = apply_mat_gene_gate([_result(P1_DOC, split_locus={"core_gene": "sexM"})], SPECS)
    assert w.withheld_reason == WITHHELD_PARALOG_CLASS


def test_a_call_without_a_paralog_class_is_judged_by_the_gate_as_before():
    from MATPredict.detect.mat_gene_gate import apply_mat_gene_gate
    doc = dict(P1_DOC, paralog_class=None)
    kept, withheld = apply_mat_gene_gate([_result(doc)], SPECS)
    assert len(kept) == 1 and withheld == []  # two flanks >= 40%


# --- end to end ----------------------------------------------------------------

def _setup_p1_locus(tmp_path):
    from tests.detect.test_idiomorph_classifier import _setup
    _setup(tmp_path)
    cdir = tmp_path / "P" / "classifiers" / "MAT"
    (cdir / "paralogs").mkdir()
    _hmm(P1, "P1", cdir / "paralogs" / "P1.hmm")
    (cdir / "manifest.yaml").write_text(yaml.safe_dump(
        {"family": "P:MAT", "paralog_classes": [{"name": "P1", "hmm": "paralogs/P1.hmm"}]}))
    cds = "".join(CODON[a] for a in P1) + "TAA"
    genome = "A" * 999 + cds + "A" * (7000 - 999 - len(cds))
    (tmp_path / "genome.fa").write_text(f">c1\n{genome}\n")
    return 1000, 999 + len(cds) - 3


def _run_p1(tmp_path):
    from MATPredict.detect.pipeline import run_pipeline
    start, end = _setup_p1_locus(tmp_path)
    spans = {"sexP": (start, end), "sexM": (start, end), "tptA": (3000, 4000), "rnhA": (5000, 6000)}
    roles = {"sexP": "core_MAT", "sexM": "core_MAT", "tptA": "flanking_conserved",
             "rnhA": "flanking_conserved"}

    def localize(*a, **k):
        return [SearchHit(KEY, g, roles[g], "c1", *spans[g], "+", i, "rec1", "tblastn_genome",
                          bitscore=b)
                for g, i, b in (("sexM", 36.1, 60.0), ("sexP", 34.0, 55.0),
                                ("tptA", 90.0, 400.0), ("rnhA", 90.0, 400.0))]

    def model(gene_name, method):
        identity = 35.0 if roles[gene_name] == "core_MAT" else 90.0
        m = _model(gene_name, "c1", *spans[gene_name], identity=identity, family_key=KEY,
                   role=roles[gene_name], method=method)
        return dataclasses.replace(m, score={"sexM": 300.0, "sexP": 150.0}.get(gene_name, 500.0))

    from MATPredict.detect import classifier as mod
    mod._CACHE.clear()
    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa", search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **kw: model(gene_name, "miniprot_refine"),
    )


def test_a_p1_locus_is_withheld_and_reported(tmp_path):
    from MATPredict.detect.mat_gene_gate import WITHHELD_PARALOG_CLASS
    out = _run_p1(tmp_path)
    assert out.results == []
    assert out.suppressed_paralog_class == 1
    (s,) = [s for s in out.suppressed_loci if s.withheld_reason == WITHHELD_PARALOG_CLASS]
    assert s.idiomorph_classifier["paralog_class"] == "P1"


# --- the build ------------------------------------------------------------------

def test_the_build_makes_paralog_hmms_from_curated_sources(tmp_path):
    from MATPredict.detect.classifier_build import build_paralog_classes
    out = tmp_path / "clf"
    (out / "paralogs").mkdir(parents=True)
    _hmm(SEXP, "sexP", out / "sexP.hmm")
    _hmm(SEXM, "sexM", out / "sexM.hmm")
    (out / "paralogs" / "P1.faa").write_text(f">GCA_000697295.1_MucIndB7402-1.0\n{P1}\n")
    (out / "paralogs" / "P1.yaml").write_text(yaml.safe_dump({
        "source": "GCA_000697295.1 (Mucor indicus B7402), BFD", "reason": "present in both",
        "evidence": "results/2026-09-29_sexM_like_paralog/NOTE.md"}))
    block = build_paralog_classes(out, {"sexP": "Plus", "sexM": "Minus"}, min_margin=25.0,
                                  check_proteins={"refP": ("sexP", SEXP), "refM": ("sexM", SEXM)})
    (entry,) = block
    assert entry["name"] == "P1" and entry["hmm"] == "paralogs/P1.hmm"
    assert (out / "paralogs" / "P1.hmm").exists()
    assert entry["n_sequences"] == 1 and entry["provenance"]["source"].startswith("GCA_000697295")
    assert entry["training_mat_proteins_classed_paralog"] == 0
    assert entry["training_mat_proteins_checked"] == 2
    assert entry["self_check_classed_paralog"] == 1


def test_a_paralog_source_without_provenance_is_refused(tmp_path):
    from MATPredict.detect.classifier_build import build_paralog_classes
    out = tmp_path / "clf"
    (out / "paralogs").mkdir(parents=True)
    _hmm(SEXP, "sexP", out / "sexP.hmm")
    _hmm(SEXM, "sexM", out / "sexM.hmm")
    (out / "paralogs" / "P1.faa").write_text(f">x\n{P1}\n")
    with pytest.raises(RuntimeError):
        build_paralog_classes(out, {"sexP": "Plus", "sexM": "Minus"}, min_margin=25.0,
                              check_proteins={})


def test_the_shipped_mucoromycota_classifier_carries_p1():
    """The real manifest lists P1, built by the script, and the shipped HMMs
    class the P1 source as P1 but no curated sexM or sexP."""
    from MATPredict.detect import classifier as mod
    from MATPredict.detect.classifier import classify, load_classifier
    from MATPredict.detect.family_registry import load_all_families
    from tests.detect.test_idiomorph_classifier import DB
    mod._CACHE.clear()
    fam = next(f for f in load_all_families(DB)
               if (f.key.phylum, f.key.locus_name) == ("Mucoromycota", "MAT"))
    clf = load_classifier(fam.idiomorph_classifier, fam)
    assert set(clf.paralogs) == {"P1"}
    assert classify(clf, [P1]).paralog_class == "P1"
    assert classify(clf, [SEXM]).paralog_class is None
    assert classify(clf, [SEXP]).paralog_class is None
    entry = next(p for p in mod.read_manifest(clf.directory)["paralog_classes"] if p["name"] == "P1")
    assert entry["training_mat_proteins_classed_paralog"] == 0
