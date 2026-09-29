"""Profile-HMM idiomorph classifier on the MODELLED proteins.

Curator's ruling, J. Stajich, 2026-09-26: add an HMM classifier as a general
step after gene modelling -- "a cleaner way to make these assignments".
Measured first (results/2026-09-26_sexMP_hmm): on modelled proteins sexM and
sexP separated 38/38 leave-one-genus-out with a worst margin of 36.6 bits
(full-length HMM) against 5.1 for blastp. The shipped Mucoromycota:MAT build
(curator option (a): curated records + UFBoot-99 sexP-clade members, Zygo 23
held out) scores 85/85 leave-one-genus-out and 23/23 on the held-out Zygo
proteins (results/2026-09-26_hmm_classifier).

Contract:
* each modelled core protein of an idiomorph-specific gene is scored against
  one HMM per idiomorph; the higher score decides `idiomorph`;
* below the family's `min_margin` the call is `undetermined`, both scores
  reported, and no hit changes;
* a decisive verdict also makes the pair's live/superseded state agree;
* no classifier, or nothing modelled: the existing decision stands.
"""
import dataclasses
from pathlib import Path

import pytest

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _model

DB = Path(__file__).resolve().parents[2] / "db"
KEY = FamilyKey("P", "MAT")
CODON = {"A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA", "E": "GAA",
         "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTT", "K": "AAA", "M": "ATG", "F": "TTT",
         "P": "CCT", "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT"}


def _curated(record, gene):
    text = (DB / "Mucoromycota" / "Mucorales" / record / "proteins.faa").read_text()
    lines = text.splitlines()
    i = next(n for n, l in enumerate(lines) if f"name={gene}|" in l)
    return lines[i + 1].strip()


SEXP = _curated("64495_cbs346-36_MAT_Plus", "sexP")
SEXM = _curated("64495_cbs110-17_MAT_Minus", "sexM")


def _hmm(seq, name, path):
    import pyhmmer
    a = pyhmmer.easel.Alphabet.amino()
    msa = pyhmmer.easel.TextMSA(name=name.encode(), sequences=[
        pyhmmer.easel.TextSequence(name=b"a", sequence=seq),
        pyhmmer.easel.TextSequence(name=b"b", sequence=seq)]).digitize(a)
    hmm, _, _ = pyhmmer.plan7.Builder(a).build_msa(msa, pyhmmer.plan7.Background(a))
    with open(path, "wb") as fh:
        hmm.write(fh)


def _setup(tmp_path, *, min_margin=10.0, classifier=True, write_hmms=True):
    (tmp_path / "P").mkdir(parents=True)
    clf = (f"    idiomorph_classifier: {{type: hmm, dir: classifiers/MAT, min_margin: {min_margin}}}\n"
           if classifier else "")
    (tmp_path / "P" / "order.yml").write_text(
        "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
        "    idiomorph_values: [Plus, Minus]\n    taxonomic_scope: [1]\n"
        "    model_idiomorph_alternatives: true\n" + clf +
        "    genes:\n"
        "      - {name: sexP, role: core_MAT, present_in_idiomorphs: [Plus]}\n"
        "      - {name: sexM, role: core_MAT, present_in_idiomorphs: [Minus]}\n"
        "      - {name: tptA, role: flanking_conserved}\n"
        "      - {name: rnhA, role: flanking_conserved}\n"
    )
    rec = tmp_path / "P" / "Fam" / "rec1"
    rec.mkdir(parents=True)
    rec.joinpath("metadata.yaml").write_text(
        "record_id: rec1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n")
    cdir = tmp_path / "P" / "classifiers" / "MAT"
    cdir.mkdir(parents=True)
    (cdir / "manifest.yaml").write_text("family: P:MAT\n")
    if write_hmms:
        _hmm(SEXP, "sexP", cdir / "sexP.hmm")
        _hmm(SEXM, "sexM", cdir / "sexM.hmm")
    # The HMG gene at 1000 encodes the SEXP protein.
    cds = "".join(CODON[a] for a in SEXP) + "TAA"
    genome = "A" * 999 + cds + "A" * (7000 - 999 - len(cds))
    (tmp_path / "genome.fa").write_text(f">c1\n{genome}\n")
    return 1000, 999 + len(cds) - 3


def _run(tmp_path, stray=False, core_models=True, flank_identity=90.0, **kw):
    start, end = _setup(tmp_path, **kw)
    spans = {"sexP": (start, end), "sexM": (start, end), "tptA": (3000, 4000),
             "rnhA": (5000, 6000)}
    roles = {"sexP": "core_MAT", "sexM": "core_MAT", "tptA": "flanking_conserved",
             "rnhA": "flanking_conserved"}

    def hit(gene, identity, bits):
        s, e = spans[gene]
        return SearchHit(KEY, gene, roles[gene], "c1", s, e, "+", identity, "rec1",
                         "tblastn_genome", bitscore=bits)

    def localize(*a, **k):
        hits = [hit("sexM", 36.1, 60.0), hit("sexP", 34.0, 55.0),
                hit("tptA", 90.0, 400.0), hit("rnhA", 90.0, 400.0)]
        if stray:  # a separate weak sexM-like hit, away from the scored HMG gene
            hits.append(SearchHit(KEY, "sexM", "core_MAT", "c1", 6500, 6600, "+", 43.0, "rec1",
                                  "tblastn_genome", bitscore=30.0))
        return hits

    def model(gene_name, method):
        if not core_models and roles[gene_name] == "core_MAT":
            return None  # no core protein models (the fragment path)
        # Flank models match their 90% hits (a real locus); the MAT-gene gate
        # (2026-09-28) counts flanks modelled at >= 40% as support.
        identity = 35.0 if roles[gene_name] == "core_MAT" else flank_identity
        m = _model(gene_name, "c1", *spans[gene_name], identity=identity, family_key=KEY,
                   role=roles[gene_name], method=method)
        # The model scores say sexM; the protein the locus encodes is sexP.
        return dataclasses.replace(m, score={"sexM": 300.0, "sexP": 150.0}.get(gene_name, 500.0))

    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **kw: model(gene_name, "miniprot_refine"),
    )


def test_without_a_classifier_the_model_scores_decide(tmp_path):
    (r,) = _run(tmp_path, classifier=False).results
    assert r.idiomorph == "Minus"
    assert r.idiomorph_classifier is None


def test_the_classifier_decides_on_the_modelled_protein(tmp_path):
    (r,) = _run(tmp_path).results
    assert r.idiomorph == "Plus"
    c = r.idiomorph_classifier
    assert c["method"] == "hmm" and c["idiomorph"] == "Plus"
    assert c["scores"]["Plus"] > c["scores"]["Minus"] + 10
    assert c["margin"] == pytest.approx(c["scores"]["Plus"] - c["scores"]["Minus"], abs=0.2)
    assert r.idiomorph_candidates[0]["idiomorph"] == "Plus"
    assert r.idiomorph_candidates[0]["basis"] == "hmm_classifier"


def test_a_decisive_verdict_makes_the_pair_agree(tmp_path):
    (r,) = _run(tmp_path).results
    assert "sexP" in r.genes_found and "sexM" not in r.genes_found
    (e,) = [e for e in r.idiomorph_resolutions if {e.winner, e.loser} == {"sexM", "sexP"}]
    assert (e.winner, e.loser, e.basis) == ("sexP", "sexM", "hmm_classifier")
    assert r.polished_genes == 3          # sexP, tptA, rnhA -- the loser never counts


def test_below_min_margin_the_call_is_undetermined_and_no_hit_changes(tmp_path):
    (r,) = _run(tmp_path, min_margin=1e6).results
    assert r.idiomorph == "undetermined"
    assert r.idiomorph_classifier["idiomorph"] == "undetermined"
    assert set(r.idiomorph_classifier["scores"]) == {"Plus", "Minus"}
    assert "sexM" in r.genes_found        # the model-pair decision's winner stands


def test_a_roster_that_names_a_missing_hmm_is_an_error(tmp_path):
    from MATPredict.detect.classifier import ClassifierError
    with pytest.raises(ClassifierError):
        _run(tmp_path, write_hmms=False)


def test_a_model_based_verdict_says_so(tmp_path):
    (r,) = _run(tmp_path).results
    assert r.idiomorph_classifier["classifier_input"] == "model"


def test_with_no_core_model_the_classifier_scores_the_hsp_fragment(tmp_path):
    """Curator's ruling 2026-09-27: when no core protein is modelled, the
    classifier scores the tblastn HSP translation instead of the identity
    comparison deciding (which called the sexP locus Minus: sexM 36.1% vs
    sexP 34.0%). Measured cases: Umbelopsis sp. M5902 (sexP +75.3 on its
    fragment, labelled Minus) and Mucor hiemalis gzMucHiem1 (sexM -46.3,
    labelled Plus)."""
    (r,) = _run(tmp_path, core_models=False).results
    assert r.idiomorph == "Plus"
    c = r.idiomorph_classifier
    assert c["classifier_input"] == "hsp_fragment"
    assert c["idiomorph"] == "Plus"
    assert c["scores"]["Plus"] > c["scores"]["Minus"] + 10


def test_a_fragment_verdict_keeps_the_min_margin(tmp_path):
    (r,) = _run(tmp_path, core_models=False, min_margin=1e6).results
    assert r.idiomorph == "undetermined"
    assert r.idiomorph_classifier["classifier_input"] == "hsp_fragment"


def test_a_fragment_call_without_flank_support_is_withheld_by_the_mat_gene_gate(tmp_path):
    """Curator's ruling 2026-09-28: absolute scores do not separate MAT genes
    from HMG paralogs on fragments, so a fragment-typed call needs >= 2 roster
    flanks modelled at >= 40% (`mat_gene_gate`)."""
    out = _run(tmp_path, core_models=False, flank_identity=35.0)
    assert out.results == []
    (s,) = [x for x in out.suppressed_loci if x.withheld_reason == "mat_gene_gate"]
    assert s.withheld_detail["classifier_input"] == "hsp_fragment"
    assert s.withheld_detail["supporting_flanks"] == []


def test_without_a_classifier_no_core_model_falls_back_to_identity(tmp_path):
    (r,) = _run(tmp_path, core_models=False, classifier=False).results
    assert r.idiomorph_classifier is None
    assert r.idiomorph == "Minus"   # the identity comparison, unchanged


def test_the_report_carries_the_classifier_verdict(tmp_path):
    from MATPredict.detect.report import _result_doc
    (r,) = _run(tmp_path).results
    doc = _result_doc(r)
    assert doc["idiomorph_classifier"]["idiomorph"] == "Plus"


def test_the_shipped_mucoromycota_classifier_is_wired_and_separates_the_references():
    """The real roster points at real files, and the built HMMs call a curated
    sexP Plus and a curated sexM Minus with the roster's min_margin."""
    from MATPredict.detect.classifier import classify, load_classifier
    from MATPredict.detect.family_registry import load_all_families
    fam = next(f for f in load_all_families(DB) if (f.key.phylum, f.key.locus_name) == ("Mucoromycota", "MAT"))
    assert fam.idiomorph_classifier and fam.idiomorph_classifier["type"] == "hmm"
    clf = load_classifier(fam.idiomorph_classifier, fam)
    assert classify(clf, [SEXP]).idiomorph == "Plus"
    assert classify(clf, [SEXM]).idiomorph == "Minus"


def test_only_mucoromycota_has_a_classifier_so_far():
    from MATPredict.detect.family_registry import load_all_families
    with_clf = {(f.key.phylum, f.key.locus_name) for f in load_all_families(DB) if f.idiomorph_classifier}
    assert with_clf == {("Mucoromycota", "MAT")}


def test_the_build_training_set_is_curated_records_plus_the_recorded_extra():
    """Curator option (a): curated records + the recorded sexP-clade extra
    file; sexM gets no extra members (its tree group is not supported)."""
    from MATPredict.detect.classifier_build import training_set
    from MATPredict.detect.family_registry import load_all_families
    fam = next(f for f in load_all_families(DB) if (f.key.phylum, f.key.locus_name) == ("Mucoromycota", "MAT"))
    rows = training_set(DB, fam, DB / "Mucoromycota" / "classifiers" / "MAT")
    by = {}
    for r in rows:
        by.setdefault((r["gene"], r["source"]), 0)
        by[(r["gene"], r["source"])] += 1
    assert by.get(("sexM", "training_extra"), 0) == 0
    assert by[("sexM", "curated_record")] >= 9
    assert by[("sexP", "training_extra")] == 70
    assert not any("zygo" in r["id"].lower() for r in rows)


def test_identity_spread_reports_unique_count_and_min_median_max(tmp_path):
    from MATPredict.detect.classifier_build import identity_spread
    aln = tmp_path / "a.afa"
    aln.write_text(">a\nMKV-LL\n>b\nMKV-LL\n>c\nMRVALL\n")
    uniq, spread = identity_spread(aln)
    assert uniq == 2
    assert spread[0] == 80.0 and spread[2] == 100.0


def test_a_verdict_changes_only_hits_at_the_scored_gene(tmp_path):
    """Found on the first 293-genome run: superseding EVERY loser hit in the
    cluster also removed separate weak hits elsewhere, dropped 11 real calls
    below the fraction floor, and they vanished from the report. Only hits
    overlapping the scored models may change."""
    (r,) = _run(tmp_path, stray=True).results
    assert r.idiomorph == "Plus"
    assert "sexM" in r.genes_found, "the stray sexM hit away from the scored gene must stay live"


def test_a_model_off_frame_by_one_base_is_still_translated(tmp_path):
    """A polished model whose exon starts one base before the codon boundary
    must still yield the real protein (found on three Mucor circinelloides
    sexM models that scored ~0 against both HMMs)."""
    from Bio.Seq import Seq
    from MATPredict.detect.pipeline import _translate_model
    cds = "".join(CODON[a] for a in SEXP)
    (tmp_path / "g.fa").write_text(f">c1\nG{cds}TAA\n")
    m = _model("sexP", "c1", 1, len(cds) + 1, family_key=KEY)   # starts one base early
    assert _translate_model(tmp_path / "g.fa", m, 1, {}) == SEXP


def test_each_gene_position_is_classified_on_its_own_protein(tmp_path):
    """Two separate HMG genes in one locus, one encoding the sexM protein and
    one the sexP protein. The pooled locus verdict must not overwrite the
    position whose own protein says otherwise (found on Radiomyces
    spectabilis and Rhizomucor pusillus: applying the pooled verdict to both
    genes removed one name and dropped the call below the fraction floor)."""
    _setup(tmp_path)
    m_cds = "".join(CODON[a] for a in SEXM) + "TAA"
    p_cds = "".join(CODON[a] for a in SEXP) + "TAA"
    a0, b0 = 1000, 1000 + len(m_cds) + 500
    genome = ("A" * (a0 - 1) + m_cds + "A" * 500 + p_cds).ljust(9000, "A")
    (tmp_path / "genome.fa").write_text(f">c1\n{genome}\n")
    spans = {"sexM": (a0, a0 + len(m_cds) - 4), "sexP": (b0, b0 + len(p_cds) - 4),
             "tptA": (7000, 7500), "rnhA": (8000, 8500)}
    roles = {"sexP": "core_MAT", "sexM": "core_MAT", "tptA": "flanking_conserved",
             "rnhA": "flanking_conserved"}

    def hit(gene, span, identity):
        return SearchHit(KEY, gene, roles[gene], "c1", *span, "+", identity, "rec1",
                         "tblastn_genome", bitscore=identity)

    def localize(*a, **k):
        # At each HMG gene both names hit; the right one wins the first pass.
        return [hit("sexM", spans["sexM"], 60.0), hit("sexP", spans["sexM"], 30.0),
                hit("sexP", spans["sexP"], 60.0), hit("sexM", spans["sexP"], 30.0),
                hit("tptA", spans["tptA"], 90.0), hit("rnhA", spans["rnhA"], 90.0)]

    def model(gene_name, method):
        return _model(gene_name, "c1", *spans[gene_name], identity=60.0, family_key=KEY,
                      role=roles[gene_name], method=method)

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **kw: model(gene_name, "miniprot_refine"),
    )
    (r,) = outcome.results
    assert {"sexM", "sexP"} <= set(r.genes_found)
