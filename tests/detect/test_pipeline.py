from __future__ import annotations
from pathlib import Path

from MATPredict.detect import pipeline as pipeline_module
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.polish import (
    STATUS_NOT_POLISH_CANDIDATE,
    STATUS_UNPOLISHED,
    ExonSpan,
    PolishModel,
)
from MATPredict.detect.search import SearchHit
from MATPredict.detect.pipeline import EvidenceFloor, run_pipeline

FAMILY = Family(FamilyKey("P", "aLocus"), "pattern", None, "^a[0-9]+$",
                [{"name": "mfa1", "role": "core_MAT"}, {"name": "pra1", "role": "core_MAT"}], [1])

ORDER_YML = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
)


def _write_order(tmp_path, text=ORDER_YML):
    (tmp_path / "P").mkdir(exist_ok=True)
    (tmp_path / "P" / "order.yml").write_text(text)


def _write_record(tmp_path, record_id="rec1", locus_name="aLocus", proteins: str = ""):
    """A curated record dir, so load_record_families can map record_id -> family."""
    record_dir = tmp_path / "P" / "Fam" / record_id
    record_dir.mkdir(parents=True, exist_ok=True)
    (record_dir / "metadata.yaml").write_text(
        f"record_id: {record_id}\n"
        f"mating_type: {{locus_name: {locus_name}, idiomorphs: [a1]}}\n"
    )
    if proteins:
        (record_dir / "proteins.faa").write_text(proteins)
    return record_dir


def _no_polish(**kwargs):
    """A polishing tool that produces no model for any gene in any window."""
    return None


def _no_localize(*args, **kwargs):
    """A localization search that finds nothing -- injected wherever a fixture
    has a routed family with zero fast-path hits, or with a partial foothold
    still missing a core_MAT gene, so neither batched rescue can reach a real
    tblastn binary."""
    return []


def _model(gene_name, contig, start, end, *, identity=80.0, strand="+",
           family_key=FAMILY.key, role="core_MAT", record_id="rec1",
           method="exonerate_refine", exons=None):
    """A single-exon PolishModel stand-in for a polishing tool's output."""
    return PolishModel(
        gene_name=gene_name, family_key=family_key, role=role, contig=contig,
        start=start, end=end, strand=strand, exons=exons or [ExonSpan(start, end)],
        identity=identity, reference_record_id=record_id, method=method,
    )


def test_run_pipeline_end_to_end_with_stubbed_search(tmp_path):
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
                SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa",
        proteome_fasta=tmp_path / "proteome.faa",
        taxid=None,
        db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
    )
    assert len(outcome.results) == 1
    result = outcome.results[0]
    assert result.family_key == FAMILY.key
    assert result.confidence == "high"
    assert result.genes_missing == []
    assert outcome.not_detected == []
    assert outcome.families_attempted == [FAMILY.key]


def _tblastn(gene_name, contig, start, end, identity=70.0, family_key=FAMILY.key):
    return SearchHit(family_key, gene_name, "core_MAT", contig, start, end, "+",
                     identity, "rec1", "tblastn_genome")


def test_genome_only_path_uses_search_localize_not_search_genomic(tmp_path):
    """The genome-only path must call the injected search_localize, and must
    never call search_genomic at all -- search_genomic's old whole-genome and
    windowed-relaxed roles are fully retired, so the name is not even bound in
    pipeline.py any more and cannot be reached from it."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    localize_calls = []

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        localize_calls.append(genome_fasta)
        return [_tblastn("mfa1", "c1", 100, 200), _tblastn("pra1", "c1", 300, 400)]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=fake_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    assert localize_calls == [tmp_path / "genome.fa"]  # exactly one batched call
    assert not hasattr(pipeline_module, "search_genomic")
    assert outcome.results[0].genes_found == ["mfa1", "pra1"]


def test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window(tmp_path):
    """A fast-path cluster missing one core gene is polished by
    polish_with_exonerate/polish_with_miniprot against a window padded around
    the EXISTING cluster's span.

    Finding 5 update: this windowed rescue is no longer the ONLY thing done for
    a partial-foothold family -- the same missing gene is also included in the
    batched genome-wide localization rescue (asserted here, and exercised
    further in
    `test_partial_foothold_familys_missing_gene_is_rescued_genome_wide`).
    Before that change, a family with PARTIAL annotation coverage got less
    search than one with none."""
    _write_order(tmp_path)
    # A 100-aa curated mfa1 -> padding = 100 * 3 * 2.0 + 2000 = 2600 bp per side.
    _write_record(
        tmp_path,
        proteins=">rec1|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 100 + "\n",
    )
    localize_calls = []
    polish_calls = []

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                          "diamond_proteome")]

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        localize_calls.append([f.key for f in families])
        return []

    def fake_exonerate(*, gene_name, window, **kwargs):
        polish_calls.append(("exonerate", gene_name, window))
        return _model("mfa1", "c1", 150, 260)

    def fake_miniprot(*, gene_name, window, **kwargs):
        polish_calls.append(("miniprot", gene_name, window))
        return _model("mfa1", "c1", 150, 260, method="miniprot_refine")

    outcome = run_pipeline(
        # Explicitly permissive: this test is about polish/segment behavior, not
        # about the admission bar, and its fixtures build single-gene clusters
        # that the curator-ruled default floor (>=2 genes) deliberately rejects.
        evidence_floor=EvidenceFloor(min_hits=1, require_core_role=False),
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=fake_localize,
        polish_with_exonerate=fake_exonerate, polish_with_miniprot=fake_miniprot,
    )
    # exactly ONE batched localization call, for the partial-foothold family
    assert localize_calls == [[FAMILY.key]]
    # only the missing gene is polished, by both tools, against the same window
    assert [c[0] for c in polish_calls] == ["exonerate", "miniprot"]
    assert {c[1] for c in polish_calls} == {"mfa1"}
    # window is the existing cluster's span (300-400) padded by 2600 either side
    assert {c[2] for c in polish_calls} == {("c1", 1, 3000)}
    result = outcome.results[0]
    assert result.genes_found == ["mfa1", "pra1"]
    # a successful polish is not "unpolished", so nothing caps the tier
    assert result.confidence == "high"
    evidence = {e.gene_name: e for e in result.gene_evidence}
    assert (evidence["mfa1"].start, evidence["mfa1"].end) == (150, 260)
    assert evidence["mfa1"].method == "exonerate_refine"
    # pra1 was found directly via the fast-path diamond hit and never
    # entered the polish stage -- distinct from an attempted-and-failed
    # STATUS_UNPOLISHED gene.
    assert evidence["pra1"].status == STATUS_NOT_POLISH_CANDIDATE


def test_gene_found_directly_never_polished_gets_not_polish_candidate_status(tmp_path):
    """A gene found directly via a confident fast-path diamond hit -- never
    localized by tblastn, never one of its family's own missing core_MAT
    genes needing rescue -- never enters the localize/polish pipeline at all.
    Its GeneEvidence.status must be STATUS_NOT_POLISH_CANDIDATE, distinct
    from STATUS_UNPOLISHED (reserved for a gene that WAS sent through both
    polish tools but that neither could confirm)."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
                SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
    )
    result = outcome.results[0]
    evidence = {e.gene_name: e for e in result.gene_evidence}
    assert evidence["mfa1"].status == STATUS_NOT_POLISH_CANDIDATE
    assert evidence["pra1"].status == STATUS_NOT_POLISH_CANDIDATE
    assert evidence["mfa1"].alternate_model is None


def test_localized_gene_neither_tool_confirms_gets_unpolished_status_not_not_polish_candidate(tmp_path):
    """A gene that WAS localized (tblastn) and sent through both polish tools,
    but neither `exonerate --refine` nor `miniprot` could produce a usable
    model, keeps the STATUS_UNPOLISHED label -- it must NOT be reported as
    STATUS_NOT_POLISH_CANDIDATE, which is reserved for a gene that never
    entered the polish stage at all."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    coords = {"mfa1": (100, 200), "pra1": (300, 400)}

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [_tblastn(name, "c1", *span) for name, span in coords.items()]

    def fake_exonerate(*, gene_name, **kwargs):
        if gene_name == "pra1":
            return None
        return _model(gene_name, "c1", *coords[gene_name])

    def fake_miniprot(*, gene_name, **kwargs):
        if gene_name == "pra1":
            return None
        return _model(gene_name, "c1", *coords[gene_name], method="miniprot_refine")

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=fake_localize,
        polish_with_exonerate=fake_exonerate, polish_with_miniprot=fake_miniprot,
    )
    evidence = {e.gene_name: e for e in outcome.results[0].gene_evidence}
    assert evidence["pra1"].status == STATUS_UNPOLISHED
    assert evidence["pra1"].status != STATUS_NOT_POLISH_CANDIDATE


def test_any_gene_unpolished_tiering_unaffected_by_not_polish_candidate_status(tmp_path):
    """`_any_gene_unpolished` (tiering) reads `PolishOutcome.status` from
    `polish_by` directly -- a case-1 gene (found directly, never a polish
    candidate) never has a `polish_by` entry at all, so it structurally
    cannot cap the tier the way a genuine `STATUS_UNPOLISHED` gene does.
    A family whose every gene is case 1 (all found directly via the
    fast-path diamond hit, nothing localized or rescued) reaches High --
    it is NOT capped to Medium the way `test_unpolished_gene_caps_tier_at_medium`
    shows a genuine STATUS_UNPOLISHED gene caps it -- proving the new
    STATUS_NOT_POLISH_CANDIDATE label has zero effect on tiering."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
                SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
    )
    result = outcome.results[0]
    evidence = {e.gene_name: e for e in result.gene_evidence}
    assert evidence["mfa1"].status == STATUS_NOT_POLISH_CANDIDATE
    assert evidence["pra1"].status == STATUS_NOT_POLISH_CANDIDATE
    # not capped -- confirms STATUS_NOT_POLISH_CANDIDATE has no tiering effect,
    # in contrast with test_unpolished_gene_caps_tier_at_medium's genuine
    # STATUS_UNPOLISHED case, which caps the identical family at Medium.
    assert result.confidence == "high"


def test_polished_agree_and_disagree_produce_identical_tier(tmp_path):
    """Two otherwise-identical scenarios -- one where the two polish tools agree
    on a gene's boundaries and one where they disagree on the SAME gene -- must
    produce the same confidence tier and the same canonical coordinates.
    Agreement is reported, never consulted by assign_tier."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    coords = {"mfa1": (100, 200), "pra1": (300, 400)}

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [_tblastn(name, "c1", *span) for name, span in coords.items()]

    def run(miniprot_shift):
        def fake_exonerate(*, gene_name, **kwargs):
            return _model(gene_name, "c1", *coords[gene_name])

        def fake_miniprot(*, gene_name, **kwargs):
            start, end = coords[gene_name]
            return _model(gene_name, "c1", start + miniprot_shift, end + miniprot_shift,
                          method="miniprot_refine")

        return run_pipeline(
            genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
            db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
            search_localize=fake_localize,
            polish_with_exonerate=fake_exonerate, polish_with_miniprot=fake_miniprot,
        )

    agree = run(0).results[0]          # both tools within the 10bp tolerance
    disagree = run(500).results[0]     # miniprot 500bp away -> polished_disagree

    assert agree.confidence == disagree.confidence == "high"
    assert [(e.gene_name, e.start, e.end) for e in agree.gene_evidence] == [
        (e.gene_name, e.start, e.end) for e in disagree.gene_evidence
    ]


def test_gene_evidence_status_and_alternate_model_reflect_polish_outcome(tmp_path):
    """`GeneEvidence.status` follows `PolishOutcome.status` per gene: a
    disagreeing gene reports the non-canonical (miniprot) model as
    `alternate_model`, an agreeing gene reports no alternate at all, and an
    unpolished gene (falls back to the raw tblastn hit) is `unpolished` with
    no alternate."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    coords = {"mfa1": (100, 200), "pra1": (300, 400)}

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [_tblastn(name, "c1", *span) for name, span in coords.items()]

    def fake_exonerate(*, gene_name, **kwargs):
        if gene_name == "pra1":
            return None  # pra1 is unpolished
        return _model(gene_name, "c1", *coords[gene_name])

    def fake_miniprot(*, gene_name, **kwargs):
        if gene_name == "pra1":
            return None
        start, end = coords[gene_name]
        # mfa1: shift miniprot's model well past tolerance -> disagree
        return _model(gene_name, "c1", start + 500, end + 500, method="miniprot_refine")

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=fake_localize,
        polish_with_exonerate=fake_exonerate, polish_with_miniprot=fake_miniprot,
    )
    evidence = {e.gene_name: e for e in outcome.results[0].gene_evidence}

    mfa1 = evidence["mfa1"]
    assert mfa1.status == "polished_disagree"
    assert mfa1.alternate_model is not None
    assert mfa1.alternate_model["method"] == "miniprot_refine"
    assert mfa1.alternate_model["start"] == 600
    assert mfa1.alternate_model["end"] == 700

    pra1 = evidence["pra1"]
    assert pra1.status == "unpolished"
    assert pra1.alternate_model is None
    assert (pra1.start, pra1.end) == (300, 400)  # raw tblastn hit


def test_unpolished_gene_caps_tier_at_medium(tmp_path):
    """A family whose genes are all localized but one gene's polish outcome is
    unpolished (neither tool produced a model) reaches at most Medium -- the
    same effect the retired second_pass_used flag had. The control run, where
    that same gene does polish, reaches High, so the cap is what makes the
    difference and not the fixture."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    coords = {"mfa1": (100, 200), "pra1": (300, 400)}

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [_tblastn(name, "c1", *span) for name, span in coords.items()]

    def run(polishable):
        def polish(*, gene_name, **kwargs):
            if gene_name not in polishable:
                return None
            return _model(gene_name, "c1", *coords[gene_name])

        return run_pipeline(
            genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
            db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
            search_localize=fake_localize,
            polish_with_exonerate=polish, polish_with_miniprot=polish,
        ).results[0]

    assert run({"mfa1", "pra1"}).confidence == "high"
    capped = run({"pra1"})  # mfa1 localized by tblastn but modelled by neither tool
    assert capped.confidence == "medium"
    assert capped.genes_found == ["mfa1", "pra1"]  # still counted as found


def test_gene_evidence_for_unpolished_gene_uses_raw_localization_hit(tmp_path):
    """An unpolished gene's GeneEvidence carries the raw tblastn hit's own
    coordinates, identity and method -- not a fabricated or missing value --
    while its polished sibling carries the polished model's."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [_tblastn("mfa1", "c1", 100, 200, identity=61.5),
                _tblastn("pra1", "c1", 300, 400, identity=72.0)]

    def polish(*, gene_name, **kwargs):
        return None if gene_name == "mfa1" else _model("pra1", "c1", 305, 395, identity=88.0)

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    evidence = {e.gene_name: e for e in outcome.results[0].gene_evidence}
    assert (evidence["mfa1"].start, evidence["mfa1"].end) == (100, 200)
    assert evidence["mfa1"].method == "tblastn_genome"
    assert evidence["mfa1"].identity == 61.5
    assert (evidence["pra1"].start, evidence["pra1"].end) == (305, 395)
    assert evidence["pra1"].method == "exonerate_refine"


def test_polished_model_attributed_to_another_family_is_rejected(tmp_path):
    """Defence in depth against this file's historic cross-family attribution
    bug class: a model that comes back keyed to a DIFFERENT family must never
    supply this family's coordinates. The gene degrades to unpolished -- its raw
    tblastn hit stands, and the tier is capped -- rather than being credited."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [_tblastn("mfa1", "c1", 100, 200), _tblastn("pra1", "c1", 300, 400)]

    def polish(*, gene_name, **kwargs):
        if gene_name == "mfa1":
            # right gene name, wrong family -- must be rejected
            return _model("mfa1", "c1", 900, 999, family_key=FamilyKey("P", "bLocus"))
        return _model("pra1", "c1", 300, 400)

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    evidence = {e.gene_name: e for e in outcome.results[0].gene_evidence}
    assert (evidence["mfa1"].start, evidence["mfa1"].end) == (100, 200)
    assert evidence["mfa1"].method == "tblastn_genome"
    assert outcome.results[0].confidence == "medium"


def test_unpolished_in_one_cluster_does_not_cap_a_different_cluster_of_the_same_family(tmp_path):
    """Same family, two independent spatial clusters (gene duplication /
    multi-allele co-occurrence is normal at MAT loci). The c1 cluster has one
    unpolished gene and is capped at medium; the c2 cluster polishes cleanly and
    must still reach high, never inherit c1's cap."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            _tblastn("mfa1", "c1", 100, 200), _tblastn("pra1", "c1", 300, 400),
            _tblastn("mfa1", "c2", 100, 200), _tblastn("pra1", "c2", 300, 400),
        ]

    def polish(*, gene_name, window, **kwargs):
        if window[0] == "c1" and gene_name == "mfa1":
            return None  # only this cluster's mfa1 fails to polish
        return _model(gene_name, window[0], 100 if gene_name == "mfa1" else 300,
                      200 if gene_name == "mfa1" else 400)

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    by_contig = {r.contig: r for r in outcome.results}
    assert len(outcome.results) == 2
    assert by_contig["c1"].confidence == "medium"
    assert by_contig["c2"].confidence == "high"


def test_family_with_zero_fast_path_hits_is_rescued_via_localization(tmp_path):
    """The blind spot this pipeline exists to close: a gene (a short pheromone
    precursor, say) that a supplied genome annotation simply does not contain
    cannot be found by searching that annotation. A family with ZERO fast-path
    hits therefore falls through to a batched tblastn localization over the
    genome, and its rescued hits are clustered and polished exactly like any
    other localized cluster -- ending in a real detected result."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    localize_calls = []
    polish_calls = []
    coords = {"mfa1": (100, 200), "pra1": (300, 400)}

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return []

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        localize_calls.append((genome_fasta, [f.key for f in families]))
        return [_tblastn(name, "c1", *span) for name, span in coords.items()]

    def polish(*, gene_name, window, **kwargs):
        polish_calls.append((gene_name, window))
        return _model(gene_name, "c1", *coords[gene_name])

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    # exactly ONE batched localization call, covering the zero-hit families
    assert localize_calls == [(tmp_path / "genome.fa", [FAMILY.key])]
    # every rescued gene is polished, by both tools
    assert sorted(c[0] for c in polish_calls) == ["mfa1", "mfa1", "pra1", "pra1"]
    assert len(outcome.results) == 1
    result = outcome.results[0]
    assert result.family_key == FAMILY.key
    assert result.genes_found == ["mfa1", "pra1"]
    assert result.confidence == "high"
    assert outcome.not_detected == []
    evidence = {e.gene_name: e for e in result.gene_evidence}
    assert (evidence["mfa1"].start, evidence["mfa1"].end) == (100, 200)
    assert evidence["mfa1"].method == "exonerate_refine"


def test_zero_hit_rescue_is_scoped_to_the_families_it_ran_for(tmp_path):
    """The batched rescue must not leak across families: only families with no
    fast-path hit are passed to search_localize, and a hit that comes back keyed
    to a family that already had a foothold is dropped rather than grafting a
    second, unrelated location onto that family."""
    _write_order(
        tmp_path,
        ORDER_YML
        + "  - locus_name: bLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^b[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: bE, role: core_MAT}\n      - {name: bW, role: core_MAT}\n",
    )
    _write_record(tmp_path)
    _write_record(tmp_path, "recB", "bLocus")
    b_key = FamilyKey("P", "bLocus")
    localize_calls = []

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        # aLocus has a foothold; bLocus has nothing at all.
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        localize_calls.append([f.key for f in families])
        return [
            SearchHit(b_key, "bE", "core_MAT", "c2", 100, 200, "+", 70.0, "recB", "tblastn_genome"),
            SearchHit(b_key, "bW", "core_MAT", "c2", 300, 400, "+", 70.0, "recB", "tblastn_genome"),
            # a stray hit keyed to the family that was NOT part of this rescue
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c9", 100, 200, "+", 70.0, "rec1", "tblastn_genome"),
        ]

    def polish(*, gene_name, window, **kwargs):
        return _model(gene_name, window[0], 100, 200, family_key=b_key)

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    assert localize_calls == [[b_key]]  # only the zero-hit family was searched for
    by_family = {}
    for result in outcome.results:
        by_family.setdefault(result.family_key, []).append(result)
    # the stray c9 aLocus hit was dropped: aLocus is reported once, on c1 only
    assert [r.contig for r in by_family[FAMILY.key]] == ["c1"]
    assert [r.contig for r in by_family[b_key]] == ["c2"]


def test_partial_foothold_familys_missing_gene_is_rescued_genome_wide(tmp_path):
    """Finding 5 regression: a family with a PARTIAL foothold must get the same
    batched genome-wide tblastn rescue a zero-hit family gets, restricted to the
    specific core_MAT genes the proteome did not place.

    This is the project's own motivating case: the pheromone receptor is
    correctly annotated but the short pheromone precursor beside it is absent
    from the same genome's own annotation -- and is not guaranteed to sit inside
    the narrow ~+-3kb window around the annotated gene. Before this fix, a
    partially-annotated family got LESS search than an unannotated one."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    localize_calls = []

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        # pra1 annotated on c1; mfa1 entirely absent from the annotation
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                          "diamond_proteome")]

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        localize_calls.append([f.key for f in families])
        return [
            # the real mfa1, far outside the narrow window around c1:300-400
            _tblastn("mfa1", "c1", 90_000, 90_300),
            # a stray hit for a gene this family already found -- must be dropped
            # rather than grafting a second, unrelated pra1 location onto it
            _tblastn("pra1", "c9", 100, 200),
        ]

    def polish(*, gene_name, window, **kwargs):
        if gene_name == "mfa1" and window[0] == "c1" and window[1] > 80_000:
            return _model("mfa1", "c1", 90_000, 90_300)
        return None

    outcome = run_pipeline(
        # Explicitly permissive: this test is about polish/segment behavior, not
        # about the admission bar, and its fixtures build single-gene clusters
        # that the curator-ruled default floor (>=2 genes) deliberately rejects.
        evidence_floor=EvidenceFloor(min_hits=1, require_core_role=False),
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    # ONE batched call, and the partial-foothold family is in it
    assert localize_calls == [[FAMILY.key]]
    # the rescued mfa1 became its own cluster and was polished there
    by_contig_start = {(r.contig, r.start): r for r in outcome.results}
    rescued = next(r for r in outcome.results if any(e.gene_name == "mfa1" for e in r.gene_evidence))
    mfa1 = next(e for e in rescued.gene_evidence if e.gene_name == "mfa1")
    assert (mfa1.contig, mfa1.start, mfa1.end) == ("c1", 90_000, 90_300)
    assert mfa1.method == "exonerate_refine"
    # the stray c9 pra1 hit, for a gene this family already had, was dropped
    assert "c9" not in {contig for contig, _ in by_contig_start}


def test_second_independent_cluster_gets_its_own_genome_wide_rescue(tmp_path):
    """Cluster-aware rescue eligibility: rescue is keyed per
    (family, gene, cluster), never per (family, gene).

    A family legitimately has more than one real, independent locus in one
    genome -- tetrapolar species with unlinked loci, and
    homothallic/heterothallic switching-cassette species. Here cluster A (c1)
    has BOTH core genes from the annotation while cluster B (c2) has only
    pra1. Asking the family-wide question "is mfa1 missing anywhere?" answers
    "no -- it is in cluster A", which is how cluster B used to be denied the
    genome-wide look and left with only the narrow ~+-2kb window around its own
    span. Cluster B's own copy of mfa1 sits at c2:10_000, far outside that
    window, so only a genome-wide rescue can reach it.

    The second half of the test is the attribution question: the rescue is
    batched and genome-wide, so it also returns an mfa1 hit sitting on top of
    cluster A, which already has mfa1. That hit must be dropped, so cluster A
    keeps its own annotated mfa1 and is neither re-polished nor widened by a
    rescue that cluster B triggered."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    localize_calls = []
    polish_windows = []

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            # cluster A on c1: complete, both core genes annotated
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
            # cluster B on c2: an independent second locus, missing mfa1
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        localize_calls.append([f.key for f in families])
        return [
            # cluster B's real mfa1: same contig as cluster B, within max_gap of
            # it so it joins that cluster, but far outside the narrow windowed
            # rescue around c2:300-400 that is all cluster B used to get.
            _tblastn("mfa1", "c2", 10_000, 10_200),
            # a batched-call by-catch hit for mfa1 landing on cluster A, which
            # already has mfa1 -- must not be grafted onto cluster A
            _tblastn("mfa1", "c1", 150, 250),
        ]

    def polish(*, gene_name, window, **kwargs):
        polish_windows.append((gene_name, window))
        if gene_name == "mfa1" and window[0] == "c2":
            return _model("mfa1", "c2", 10_000, 10_200)
        return None

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    # The genome-wide rescue ran at all -- under family-wide eligibility mfa1
    # counted as "found" and search_localize would never have been called.
    assert localize_calls == [[FAMILY.key]]

    by_contig = {}
    for result in outcome.results:
        by_contig.setdefault(result.contig, []).append(result)

    # cluster B got mfa1 from the genome-wide rescue, in its own cluster
    assert len(by_contig["c2"]) == 1
    cluster_b = by_contig["c2"][0]
    assert sorted(cluster_b.genes_found) == ["mfa1", "pra1"]
    b_mfa1 = next(e for e in cluster_b.gene_evidence if e.gene_name == "mfa1")
    assert (b_mfa1.contig, b_mfa1.start, b_mfa1.end) == ("c2", 10_000, 10_200)
    assert b_mfa1.method == "exonerate_refine"

    # cluster A is untouched: its mfa1 is still its own annotated one, and the
    # rescue's c1 by-catch was never folded in or polished.
    assert len(by_contig["c1"]) == 1
    cluster_a = by_contig["c1"][0]
    a_mfa1 = next(e for e in cluster_a.gene_evidence if e.gene_name == "mfa1")
    assert (a_mfa1.start, a_mfa1.end) == (100, 200)
    assert a_mfa1.method == "diamond_proteome"
    assert cluster_a.confidence == "high"
    assert [w for _gene, w in polish_windows if w[0] == "c1"] == []
    # every polish window opened was for cluster B's missing gene
    assert {gene for gene, _w in polish_windows} == {"mfa1"}


#: A three-gene family, needed to build a rescue CHAIN: one gene supplies the
#: legitimate gap a rescue hit fills, a second supplies the spurious hit that
#: chains in behind it, and a third anchors the cluster.
THREE_GENE_YML = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
    "      - {name: x1, role: core_MAT}\n"
)


def test_chained_in_rescue_hit_cannot_affect_a_gene_the_cluster_already_has(tmp_path):
    """`_RescueScope.accepts` judges a rescue hit against a cluster's ORIGINAL
    span, but `cluster_hits` merges by CHAINING -- so a hit `accepts` correctly
    rules "not on top of cluster A" can still be pulled INTO cluster A by a
    second, separately-accepted rescue hit that bridges the gap.

    Layout (max_gap = 25_000):

    * cluster A on c1 has mfa1 (100-200) and x1 (300-400) from the annotation,
      and is missing pra1;
    * cluster B on c2 has pra1 (300-400) and is missing mfa1 and x1, which is
      what puts x1 into the family's rescue scope at all even though A has it;
    * the batched rescue returns a legitimate pra1 for A at c1:20_000 (accepted
      -- A genuinely lacks pra1) and a spurious x1 at c1:44_000 (also accepted:
      it is 43_600 bp past A's original end of 400, far outside the +-25_000
      guard `accepts` applies).

    The definitive clustering then chains 400 -> 20_000 (19_600) -> 44_000
    (23_900) and the spurious x1 lands inside cluster A. Before the containment
    fix that made x1 a polish candidate in A; the polish fails (it is spurious)
    and `_any_gene_unpolished` dropped A from high to medium -- entirely because
    cluster B had a gap in a DIFFERENT gene.

    After the fix, a gene with a non-localized hit already in that specific
    cluster is never a polish candidate there, so A keeps its annotated x1
    coordinates and its high tier.
    """
    _write_order(tmp_path, THREE_GENE_YML)
    _write_record(tmp_path)
    polish_windows = []

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            # cluster A on c1: has mfa1 and x1, missing pra1
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "x1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
            # cluster B on c2: an independent locus with only pra1
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            # legitimate: cluster A really is missing pra1
            _tblastn("pra1", "c1", 20_000, 20_100),
            # spurious/paralogous x1, outside A's original +-max_gap guard but
            # chained into A by the pra1 hit above
            _tblastn("x1", "c1", 44_000, 44_100),
        ]

    def polish(*, gene_name, window, **kwargs):
        polish_windows.append((gene_name, window))
        if gene_name == "pra1" and window[0] == "c1":
            return _model("pra1", "c1", 20_000, 20_100)
        return None  # the spurious x1 is modelled by neither tool

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )

    cluster_a = next(r for r in outcome.results if r.contig == "c1")
    # A's tier is not capped by a failed polish of a gene it already had
    assert sorted(cluster_a.genes_found) == ["mfa1", "pra1", "x1"]
    assert cluster_a.confidence == "high"
    # A's x1 evidence is still its own annotated diamond hit, not the spurious HSP
    a_x1 = next(e for e in cluster_a.gene_evidence if e.gene_name == "x1")
    assert (a_x1.contig, a_x1.start, a_x1.end) == ("c1", 300, 400)
    assert a_x1.method == "diamond_proteome"
    assert a_x1.status == STATUS_NOT_POLISH_CANDIDATE
    # the chained-in spurious x1 never became a polish candidate in cluster A;
    # the only gene polished on c1 is pra1, which A genuinely lacked
    assert {gene for gene, w in polish_windows if w[0] == "c1"} == {"pra1"}
    # the legitimate rescue for the gene A really was missing still happened
    a_pra1 = next(e for e in cluster_a.gene_evidence if e.gene_name == "pra1")
    assert (a_pra1.start, a_pra1.end) == (20_000, 20_100)
    assert a_pra1.method == "exonerate_refine"


def test_rescue_eligibility_is_scoped_per_family_gene_and_cluster(tmp_path):
    """Unit-level check of the eligibility rule itself, including the zero-hit
    case this change must not regress."""
    b_family = Family(FamilyKey("P", "bLocus"), "pattern", None, "^b[0-9]+$",
                      [{"name": "mfa1", "role": "core_MAT"}], [1])
    complete = SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome")
    complete2 = SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")
    partial = SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec1", "diamond_proteome")

    from MATPredict.detect.clustering import cluster_hits

    # bLocus has no hits at all -> whole gene set in scope (the original
    # zero-hit behaviour, unchanged).
    scope = pipeline_module._localization_rescue_targets(
        cluster_hits([complete, complete2, partial], max_gap=25_000), [FAMILY, b_family]
    )
    assert scope.genes_by_family[b_family.key] is None
    # aLocus is eligible only for mfa1, and only because cluster B lacks it
    assert scope.genes_by_family[FAMILY.key] == {"mfa1"}
    assert [f.key for f in scope.families] == [FAMILY.key, b_family.key]

    # a gene neither cluster is missing is never in scope
    assert not scope.accepts(_tblastn("pra1", "c9", 100, 200), 25_000)
    # cluster B's missing gene is in scope anywhere it is not already covered
    assert scope.accepts(_tblastn("mfa1", "c2", 10_000, 10_200), 25_000)
    # ...but not on top of cluster A, which already has mfa1
    assert not scope.accepts(_tblastn("mfa1", "c1", 150, 250), 25_000)
    # a far-away second mfa1 location on cluster A's contig is still allowed:
    # it cannot be folded into cluster A, so it is a candidate new locus
    assert scope.accepts(_tblastn("mfa1", "c1", 500_000, 500_200), 25_000)
    # the zero-hit family's own gene names are accepted for its own key only
    assert scope.accepts(_tblastn("mfa1", "c1", 150, 250, family_key=b_family.key), 25_000)
    # a family the rescue never ran for gains nothing
    assert not scope.accepts(_tblastn("mfa1", "c1", 150, 250, family_key=FamilyKey("P", "zLocus")), 25_000)


def test_rescued_cluster_with_an_unpolished_gene_is_capped_at_medium(tmp_path):
    """Polish eligibility is per-cluster, not gated on a global "genome-only"
    flag: a tblastn-rescued cluster on the FAST path is polished like any other
    localized cluster, so a gene neither tool can model there still caps the
    family's tier -- which a global flag would have silently prevented."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    coords = {"mfa1": (100, 200), "pra1": (300, 400)}

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return []

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [_tblastn(name, "c1", *span) for name, span in coords.items()]

    def run(polishable):
        def polish(*, gene_name, **kwargs):
            if gene_name not in polishable:
                return None
            return _model(gene_name, "c1", *coords[gene_name])

        return run_pipeline(
            genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa",
            taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
            search_fast_path=fake_fast_path, search_localize=fake_localize,
            polish_with_exonerate=polish, polish_with_miniprot=polish,
        ).results[0]

    assert run({"mfa1", "pra1"}).confidence == "high"
    capped = run({"pra1"})  # mfa1 rescued by tblastn but modelled by neither tool
    assert capped.confidence == "medium"
    assert capped.genes_found == ["mfa1", "pra1"]


def test_polished_model_on_another_contig_is_rejected(tmp_path):
    """`_own_model`'s defensive guard checks the contig too: a model whose
    coordinates name a contig the polish window was not opened on cannot belong
    to this cluster, and must degrade the gene to unpolished rather than write
    another part of the genome's coordinates into this cluster's evidence."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [_tblastn("mfa1", "c1", 100, 200), _tblastn("pra1", "c1", 300, 400)]

    def polish(*, gene_name, **kwargs):
        if gene_name == "mfa1":
            return _model("mfa1", "c7", 900, 999)  # right family and gene, wrong contig
        return _model("pra1", "c1", 300, 400)

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    evidence = {e.gene_name: e for e in outcome.results[0].gene_evidence}
    assert (evidence["mfa1"].start, evidence["mfa1"].end) == (100, 200)
    assert evidence["mfa1"].method == "tblastn_genome"
    assert outcome.results[0].confidence == "medium"


def test_segment_span_covers_every_gene_it_reports(tmp_path):
    """A rescued/polished gene can sit outside the raw cluster span frozen at
    clustering time. The reported segment must cover the union of that span and
    every gene evidence coordinate on its contig -- otherwise the GFF3 carries a
    gene feature outside its own parent MAT_locus feature's declared range."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        # only pra1 is annotated; the cluster span is frozen at 300-400
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                          "diamond_proteome")]

    def polish(*, gene_name, **kwargs):
        return _model("mfa1", "c1", 150, 260)  # rescued gene, entirely left of 300

    outcome = run_pipeline(
        # Explicitly permissive: this test is about polish/segment behavior, not
        # about the admission bar, and its fixtures build single-gene clusters
        # that the curator-ruled default floor (>=2 genes) deliberately rejects.
        evidence_floor=EvidenceFloor(min_hits=1, require_core_role=False),
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=_no_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )
    result = outcome.results[0]
    assert [(s.contig, s.start, s.end) for s in result.segments] == [("c1", 150, 400)]
    assert (result.start, result.end) == (150, 400)
    # every reported gene really is inside its own segment
    for gene in result.gene_evidence:
        segment = next(s for s in result.segments if s.contig == gene.contig)
        assert segment.start <= gene.start and gene.end <= segment.end


def test_short_orf_gene_reported_as_not_searchable_not_missing(tmp_path):
    _write_order(tmp_path)
    _write_record(
        tmp_path,
        proteins=(
            ">rec1|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 41 + "\n"
            ">rec1|gene_index=1|name=pra1|role=core_MAT\n" + "M" * 300 + "\n"
        ),
    )

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        search_localize=_no_localize,
        # mfa1 is genuinely not found, even by the polish rescue
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    assert outcome.results[0].genes_missing == []
    assert outcome.results[0].genes_not_searchable == ["mfa1"]


def test_short_orf_split_discriminates_three_buckets(tmp_path):
    """Three core genes, three distinct fates: mfa1 is short (41 aa) and
    genuinely never found -> genes_not_searchable. pra2 is normal-length and
    genuinely never found -> genes_missing. pra1 is found -> genes_found.
    A test with only one gene per bucket can't tell "correctly bucketed" from
    "coincidentally correct because nothing else is missing" -- this proves
    the split logic actually discriminates across all three cases at once."""
    three_gene_family = Family(
        FamilyKey("P", "aLocus"), "pattern", None, "^a[0-9]+$",
        [{"name": "mfa1", "role": "core_MAT"},
         {"name": "pra1", "role": "core_MAT"},
         {"name": "pra2", "role": "core_MAT"}],
        [1],
    )
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
        "      - {name: pra2, role: core_MAT}\n",
    )
    _write_record(
        tmp_path,
        proteins=(
            ">rec1|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 41 + "\n"
            ">rec1|gene_index=1|name=pra1|role=core_MAT\n" + "M" * 300 + "\n"
            ">rec1|gene_index=2|name=pra2|role=core_MAT\n" + "M" * 300 + "\n"
        ),
    )

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [SearchHit(three_gene_family.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                           "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        search_localize=_no_localize,
        # neither mfa1 nor pra2 is found, even by the polish rescue
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        ambiguity_floor=0.3,  # only 1 of 3 genes is found by design -- lower the floor
        # so this single-family cluster still clears scoring and isn't dropped,
        # without needing a fourth gene just to satisfy an unrelated threshold.
    )
    assert outcome.results[0].genes_found == ["pra1"]
    assert outcome.results[0].genes_missing == ["pra2"]
    assert outcome.results[0].genes_not_searchable == ["mfa1"]


def test_short_orf_scan_is_scoped_per_family_and_uses_the_longest_curated_protein(tmp_path):
    """Finding 2 regression, modelled on the real MATsc `cha1` case: one
    curated record holds a 22-aa FRAGMENT of a gene that is normal length
    elsewhere. Scanning globally and keeping the SHORTEST length made that
    fragment condemn the gene everywhere, falsely reporting a perfectly
    searchable gene as "not searchable by this method"."""
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n"
        "  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: cha1, role: flanking_variable}\n"
        "  - locus_name: bLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^b[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: bE, role: core_MAT}\n      - {name: cha1, role: flanking_variable}\n",
    )
    # aLocus: one record has a full-length cha1, another has only a 22-aa fragment.
    _write_record(tmp_path, "recFull", "aLocus", proteins=(
        ">recFull|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 300 + "\n"
        ">recFull|gene_index=1|name=cha1|role=flanking_variable\n" + "M" * 400 + "\n"
    ))
    _write_record(tmp_path, "recFrag", "aLocus", proteins=(
        ">recFrag|gene_index=0|name=cha1|role=flanking_variable\n" + "M" * 22 + "\n"
    ))
    # bLocus's own cha1 evidence is genuinely short everywhere.
    _write_record(tmp_path, "recB", "bLocus", proteins=(
        ">recB|gene_index=0|name=bE|role=core_MAT\n" + "M" * 300 + "\n"
        ">recB|gene_index=1|name=cha1|role=flanking_variable\n" + "M" * 20 + "\n"
    ))

    a_key, b_key = FamilyKey("P", "aLocus"), FamilyKey("P", "bLocus")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(a_key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "recFull", "diamond_proteome"),
            SearchHit(b_key, "bE", "core_MAT", "c9", 100, 200, "+", 95.0, "recB", "diamond_proteome"),
        ]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    by_family = {r.family_key: r for r in outcome.results}
    # aLocus has a full-length curated cha1 -> genuinely missing, not "unsearchable"
    assert by_family[a_key].genes_missing == ["cha1"]
    assert by_family[a_key].genes_not_searchable == []
    # bLocus's best curated cha1 really is short -> not searchable by this method
    assert by_family[b_key].genes_missing == []
    assert by_family[b_key].genes_not_searchable == ["cha1"]


def test_sub_floor_families_are_reported_as_not_detected_not_dropped(tmp_path):
    """Finding 6 regression: the spec requires a sub-floor result to be
    reported as "not detected" listing which families were attempted and why
    each fell short -- never silently omitted."""
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n"
        "  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: g1, role: core_MAT}\n      - {name: g2, role: core_MAT}\n"
        "      - {name: g3, role: core_MAT}\n      - {name: g4, role: core_MAT}\n"
        "  - locus_name: bLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^b[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: bE, role: core_MAT}\n      - {name: bW, role: core_MAT}\n",
    )
    _write_record(tmp_path, "recA", "aLocus")
    _write_record(tmp_path, "recB", "bLocus")
    a_key, b_key = FamilyKey("P", "aLocus"), FamilyKey("P", "bLocus")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        # 1 of aLocus's 4 genes -> fraction 0.25, below the 0.5 floor.
        return [SearchHit(a_key, "g1", "core_MAT", "c1", 100, 200, "+", 90.0, "recA", "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=_no_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    assert outcome.results == []
    not_detected = {n.family_key: n for n in outcome.not_detected}
    assert set(not_detected) == {a_key, b_key}
    assert not_detected[a_key].best_fraction_found == 0.25
    assert "below the ambiguity floor" in not_detected[a_key].reason
    assert not_detected[a_key].genes_found == ["g1"]
    assert not_detected[b_key].best_fraction_found == 0.0
    assert "no reference-protein hits" in not_detected[b_key].reason
    assert sorted(outcome.families_attempted, key=lambda k: k.locus_name) == [a_key, b_key]


def test_isolated_single_hit_is_low_tier(tmp_path):
    """Finding 6, second half: "low" was unreachable because score_cluster
    never emits a score for a family with zero hits. Per the spec's own Low
    description, a single isolated hit with nothing else from the family
    nearby is Low."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=_no_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    assert outcome.results[0].confidence == "low"


def test_genes_split_across_contigs_are_reported_separately_by_default(tmp_path):
    """Curator's ruling, 2026-09-20: cross-contig fragment merging is OFF by
    default. `_fragmented_family_segments`' merge can name a locus by
    whichever contig its coordinates default to -- Basidiomycota order
    testing found a real Cryptococcus deneoformans MAT locus (both genes
    matching at 100% identity) reported cross-contig as `fragmented`, low
    confidence, `idiomorph: undetermined`, with its top-level contig/start/end
    naming a chromosome that held almost none of the call's own evidence.

    The identical two-contig fixture from the merge test (mfa1 on c1, pra1 on
    c2, no flank) now reports as TWO independent single-contig calls instead
    of one merged multi-segment one -- each partial, each fragmented=False,
    each on its own real contig with its own real coordinates. Nothing is
    silently lost: a genuinely fragmented locus is still visible, just as two
    honest partial calls rather than one call that might misname itself."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\n" + "A" * 1000 + "\n>c2\n" + "A" * 1000 + "\n")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    outcome = run_pipeline(
        genome_fasta=genome, proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        search_localize=_no_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    assert len(outcome.results) == 2
    assert all(r.fragmented is False for r in outcome.results)
    by_contig = {r.contig: r for r in outcome.results}
    assert by_contig["c1"].genes_found == ["mfa1"]
    assert by_contig["c1"].start == 100 and by_contig["c1"].end == 200
    assert by_contig["c2"].genes_found == ["pra1"]
    assert by_contig["c2"].start == 300 and by_contig["c2"].end == 400


def test_cross_contig_merge_is_still_available_when_explicitly_enabled(tmp_path):
    """Finding 5: a family whose core genes land on different contigs, with no
    single cluster carrying them all, is one multi-segment locus with
    fragmented=True and a one-tier confidence downgrade -- but ONLY when the
    caller opts in with allow_cross_contig_fragments=True. The merge machinery
    itself is unchanged and still correct; only the default flipped."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\n" + "A" * 1000 + "\n>c2\n" + "A" * 1000 + "\n")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    outcome = run_pipeline(
        genome_fasta=genome, proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        # Each of the two clusters is missing one of the family's core genes,
        # so with cluster-aware rescue eligibility both are now genome-wide
        # rescue targets and the batched localization runs (finding nothing).
        search_localize=_no_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        allow_cross_contig_fragments=True,
    )
    assert len(outcome.results) == 1
    result = outcome.results[0]
    assert result.fragmented is True
    assert sorted(result.genes_found) == ["mfa1", "pra1"]
    assert [(s.contig, s.start, s.end) for s in result.segments] == [("c1", 100, 200), ("c2", 300, 400)]
    # contig_edge_distance is populated from the real genome FASTA
    assert result.segments[0].contig_edge_distance == 99
    assert result.segments[1].contig_edge_distance == 299
    # would have been "high" on one contig; multi-segment downgrades one tier
    assert result.confidence == "medium"


def test_same_family_on_two_contigs_each_complete_is_not_fragmented(tmp_path):
    """Gene duplication / multi-allele co-occurrence is normal at MAT loci.
    Two COMPLETE copies of a family on different contigs are two real loci,
    never one fragmented locus."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c2", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    assert len(outcome.results) == 2
    assert all(r.fragmented is False for r in outcome.results)


def test_fragmented_family_with_a_separate_independent_cluster_reports_both(tmp_path):
    """Finding C regression: once a family is judged fragmented (genes split
    across c1/c2, merged into one multi-segment call), a genuine SEPARATE
    above-floor cluster for that same family on a third contig (e.g. real
    duplication) must still be reported on its own -- not silently dropped
    just because it shares a family_key with the fragmented call.

    Three core genes so the merge (c1 has mfa1+pra1, c2 has pra2) covers the
    whole family without any single cluster being complete on its own -- c3
    independently has mfa1+pra1 too (2 of 3, above the ambiguity floor), a
    real second, unrelated cluster that must not be suppressed."""
    three_gene_family = Family(
        FamilyKey("P", "aLocus"), "pattern", None, "^a[0-9]+$",
        [{"name": "mfa1", "role": "core_MAT"},
         {"name": "pra1", "role": "core_MAT"},
         {"name": "pra2", "role": "core_MAT"}],
        [1],
    )
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
        "      - {name: pra2, role: core_MAT}\n",
    )
    _write_record(tmp_path)
    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\n" + "A" * 1000 + "\n>c2\n" + "A" * 1000 + "\n>c3\n" + "A" * 1000 + "\n")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            # c1 + c2 together make one fragmented call.
            SearchHit(three_gene_family.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(three_gene_family.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(three_gene_family.key, "pra2", "core_MAT", "c2", 500, 600, "+", 95.0, "rec1", "diamond_proteome"),
            # c3 independently has 2 of 3 core genes -- a real, separate, incomplete-but-above-floor locus.
            SearchHit(three_gene_family.key, "mfa1", "core_MAT", "c3", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(three_gene_family.key, "pra1", "core_MAT", "c3", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    outcome = run_pipeline(
        genome_fasta=genome, proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        # No single cluster here carries all three core genes, so every cluster
        # is a cluster-aware rescue target and the batched localization runs.
        search_localize=_no_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        ambiguity_floor=0.5,
        allow_cross_contig_fragments=True,
    )
    assert len(outcome.results) == 2
    by_fragmented = {r.fragmented: r for r in outcome.results}
    assert True in by_fragmented and False in by_fragmented
    fragmented = by_fragmented[True]
    assert sorted({s.contig for s in fragmented.segments}) == ["c1", "c2"]
    independent = by_fragmented[False]
    assert independent.contig == "c3"
    assert sorted(independent.genes_found) == ["mfa1", "pra1"]


def test_gene_evidence_prefers_a_polished_model_over_a_higher_identity_raw_hit(tmp_path):
    """Finding 4 regression: `_gene_evidence` must never rank a raw tblastn
    `pident` against an exonerate/miniprot identity as if they were the same
    number -- the spec (Stage 2) says those are not comparable across tools. A
    polished model wins over the raw hit it refines whatever the identity
    figures say.

    c1's mfa1 carries a raw tblastn hit with a high pident (99.0) and is also
    polished, with a much LOWER identity figure (60.0). Comparing the two
    numbers picks the unrefined raw hit; the correct answer is the polished
    model. The contest is settled inside c1, using only c1's own hits.

    c2 independently holds its own raw mfa1. Since Finding 5 that is a separate
    segment's separate evidence and is reported in its own right, never merged
    into or compared against c1's -- asserted here so this test pins the
    per-segment attribution too."""
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
        "      - {name: pra2, role: core_MAT}\n",
    )
    _write_record(tmp_path)
    key = FamilyKey("P", "aLocus")
    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\n" + "A" * 3000 + "\n>c2\n" + "A" * 3000 + "\n")

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            # a high raw pident for the gene that c1 also polishes, below
            SearchHit(key, "mfa1", "core_MAT", "c1", 100, 200, "+", 99.0, "rec1", "tblastn_genome"),
            SearchHit(key, "pra1", "core_MAT", "c1", 300, 400, "+", 70.0, "rec1", "tblastn_genome"),
            # the SAME gene again on the other contig -- a separate segment's own hit
            SearchHit(key, "mfa1", "core_MAT", "c2", 100, 200, "+", 99.0, "rec1", "tblastn_genome"),
            SearchHit(key, "pra2", "core_MAT", "c2", 500, 600, "+", 70.0, "rec1", "tblastn_genome"),
        ]

    def polish(*, gene_name, window, **kwargs):
        # only mfa1 in c1's window is modelled, and with a LOW identity figure
        if gene_name == "mfa1" and window[0] == "c1":
            return _model("mfa1", "c1", 120, 190, identity=60.0, family_key=key)
        return None

    outcome = run_pipeline(
        genome_fasta=genome, proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=fake_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
        allow_cross_contig_fragments=True,
    )
    fragmented = [r for r in outcome.results if r.fragmented]
    assert len(fragmented) == 1
    evidence = {(e.contig, e.gene_name): e for e in fragmented[0].gene_evidence}
    c1_mfa1 = evidence[("c1", "mfa1")]
    assert c1_mfa1.method == "exonerate_refine"
    assert (c1_mfa1.start, c1_mfa1.end) == (120, 190)
    assert c1_mfa1.identity == 60.0
    # c2's own raw mfa1 stands on its own segment, unaffected by c1's model
    c2_mfa1 = evidence[("c2", "mfa1")]
    assert c2_mfa1.method == "tblastn_genome"
    assert (c2_mfa1.start, c2_mfa1.end) == (100, 200)


def test_contig_edge_distance_populated_regardless_of_other_families_fragmentation(tmp_path):
    """Finding D regression: contig_edge_distance must be populated (or left
    None) consistently per-result, never depending on whether some OTHER
    family in the same run happened to be fragmented. Here NO family is
    fragmented, but the genome FASTA is readable, so the single-contig
    result's segment must still get a real contig_edge_distance -- not None,
    which is what the old "only read the genome when some family is
    fragmented" gate produced."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\n" + "A" * 1000 + "\n")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    outcome = run_pipeline(
        genome_fasta=genome, proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    assert len(outcome.results) == 1
    assert outcome.results[0].fragmented is False
    assert outcome.results[0].segments[0].contig_edge_distance == 99


def test_detection_result_carries_per_gene_evidence(tmp_path):
    """Finding 7: identity, coverage, coordinates, role and the matched
    curated record must survive to the DetectionResult, not be discarded."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 91.5, "rec1",
                      "diamond_proteome", coverage=77.5),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "-", 88.0, "rec1",
                      "diamond_proteome", coverage=99.0),
        ]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    evidence = {e.gene_name: e for e in outcome.results[0].gene_evidence}
    assert evidence["mfa1"].identity == 91.5
    assert evidence["mfa1"].coverage == 77.5
    assert evidence["mfa1"].reference_record_id == "rec1"
    assert evidence["pra1"].strand == "-"
    assert evidence["pra1"].start == 300
    assert outcome.results[0].reference_records == ["rec1"]


def test_pipeline_output_feeds_the_report_writers_directly(tmp_path):
    """End-to-end boundary check: whatever run_pipeline returns must be exactly
    what report.py consumes, so "not detected" entries and per-gene evidence
    actually reach the CLI's files rather than stopping at an internal type."""
    import yaml

    from MATPredict.detect.report import write_detection_gff3, write_detection_report

    _write_order(
        tmp_path,
        ORDER_YML
        + "  - locus_name: bLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^b[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: bE, role: core_MAT}\n      - {name: bW, role: core_MAT}\n",
    )
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1",
                      "diamond_proteome", coverage=80.0),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 93.0, "rec1",
                      "diamond_proteome", coverage=90.0),
        ]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=_no_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    write_detection_gff3(outcome, tmp_path / "out.gff3")
    write_detection_report(outcome, tmp_path / "out.yaml")

    gff3 = (tmp_path / "out.gff3").read_text()
    assert "Name=mfa1" in gff3 and "identity=95.0" in gff3
    doc = yaml.safe_load((tmp_path / "out.yaml").read_text())
    assert doc["detected"][0]["gene_evidence"][0]["coverage"] == 80.0
    # the family that never cleared the floor is reported, not dropped
    assert [n["family"] for n in doc["not_detected"]] == ["P:bLocus"]


def test_fragmented_segments_each_keep_their_own_evidence_for_a_shared_gene_name(tmp_path):
    """Finding 5 regression: when two segments of ONE fragmented call both
    genuinely hold a hit for the SAME gene name, neither segment's evidence may
    be discarded by a cross-cluster "best per gene name" selection.

    The fragmentation cover is a CONTRIBUTION test, not a disjointness test: a
    cluster is admitted for contributing >=1 not-yet-covered core gene, and
    nothing rejects it for also sharing a gene name with an already-chosen
    cluster. Gene duplication / multi-allele co-occurrence is normal at MAT
    loci, so two pieces of a split locus sharing a gene name is expected.

    Here c1 carries a real annotated diamond mfa1 at 100-200 plus pra1, and c2
    carries pra2 plus -- via its own legitimate windowed rescue for the core
    gene it lacks -- a polished mfa1 at 700-800. Both are real evidence for two
    different genomic copies in two different segments. Collapsing by bare gene
    name made the polished c2 model outrank c1's annotated hit unconditionally
    and drop c1's real coordinates from the report entirely, leaving segment c1
    declared 100-400 with no gene explaining its first 200 bp."""
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
        "      - {name: pra2, role: core_MAT}\n",
    )
    _write_record(tmp_path)
    key = FamilyKey("P", "aLocus")
    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\n" + "A" * 5000 + "\n>c2\n" + "A" * 5000 + "\n")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(key, "pra2", "core_MAT", "c2", 500, 600, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    def polish(*, gene_name, window, **kwargs):
        # c2's windowed rescue legitimately models the mfa1 that c2 lacks.
        if gene_name == "mfa1" and window[0] == "c2":
            return _model("mfa1", "c2", 700, 800, identity=55.0, family_key=key)
        return None

    outcome = run_pipeline(
        # Explicitly permissive: this test is about polish/segment behavior, not
        # about the admission bar, and its fixtures build single-gene clusters
        # that the curator-ruled default floor (>=2 genes) deliberately rejects.
        evidence_floor=EvidenceFloor(min_hits=1, require_core_role=False),
        genome_fasta=genome, proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=_no_localize,
        polish_with_exonerate=polish, polish_with_miniprot=polish,
        allow_cross_contig_fragments=True,
    )
    fragmented = [r for r in outcome.results if r.fragmented]
    assert len(fragmented) == 1
    result = fragmented[0]
    assert sorted(result.genes_found) == ["mfa1", "pra1", "pra2"]

    by_place = {(e.contig, e.gene_name): e for e in result.gene_evidence}
    # BOTH segments' mfa1 survive, each attributed to its own segment.
    assert ("c1", "mfa1") in by_place, "c1's real annotated mfa1 was dropped"
    assert ("c2", "mfa1") in by_place, "c2's rescued mfa1 was dropped"
    c1_mfa1 = by_place[("c1", "mfa1")]
    assert (c1_mfa1.start, c1_mfa1.end) == (100, 200)
    assert c1_mfa1.method == "diamond_proteome"
    assert c1_mfa1.status == STATUS_NOT_POLISH_CANDIDATE
    c2_mfa1 = by_place[("c2", "mfa1")]
    assert (c2_mfa1.start, c2_mfa1.end) == (700, 800)
    assert c2_mfa1.method == "exonerate_refine"
    # the single-copy genes are untouched and still attributed to their own segment
    assert (by_place[("c1", "pra1")].start, by_place[("c1", "pra1")].end) == (300, 400)
    assert (by_place[("c2", "pra2")].start, by_place[("c2", "pra2")].end) == (500, 600)
    assert len(result.gene_evidence) == 4


def test_fragmented_segments_keep_both_raw_hits_for_a_shared_gene_name(tmp_path):
    """Finding 5 regression, variant B: the same collapse with NO polishing at
    all. Two raw `diamond_proteome` mfa1 hits, one per segment, used to be
    resolved by the identity tie-break, silently dropping the lower-identity
    one even though it is the only evidence on its own segment."""
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
        "      - {name: pra2, role: core_MAT}\n",
    )
    _write_record(tmp_path)
    key = FamilyKey("P", "aLocus")
    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\n" + "A" * 5000 + "\n>c2\n" + "A" * 5000 + "\n")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(key, "mfa1", "core_MAT", "c2", 100, 200, "+", 70.0, "rec1", "diamond_proteome"),
            SearchHit(key, "pra2", "core_MAT", "c2", 500, 600, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    outcome = run_pipeline(
        genome_fasta=genome, proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=_no_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        allow_cross_contig_fragments=True,
    )
    fragmented = [r for r in outcome.results if r.fragmented]
    assert len(fragmented) == 1
    by_place = {(e.contig, e.gene_name): e for e in fragmented[0].gene_evidence}
    assert by_place[("c1", "mfa1")].identity == 95.0
    assert by_place[("c2", "mfa1")].identity == 70.0
    assert len(fragmented[0].gene_evidence) == 4


# --- per-locus cluster gap derivation ---


def _gap_recording_cluster_hits(recorded):
    """Wrap the real `cluster_hits`, recording the max_gap each call was given."""
    real = pipeline_module.cluster_hits

    def wrapper(hits, max_gap):
        recorded.append(max_gap)
        return real(hits, max_gap=max_gap)

    return wrapper


WIDE_GAP_ORDER_YML = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    max_cluster_gap_bp: 50000\n"
    "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
)


def _run_recording_gap(tmp_path, monkeypatch, order_text, **kwargs):
    _write_order(tmp_path, order_text)
    _write_record(tmp_path)
    recorded: list[int] = []
    monkeypatch.setattr(
        pipeline_module, "cluster_hits", _gap_recording_cluster_hits(recorded)
    )

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [_tblastn("mfa1", "c1", 100, 200), _tblastn("pra1", "c1", 300, 400)]

    run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None,
        taxid=kwargs.pop("taxid", None),
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=fake_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        **kwargs,
    )
    return recorded


def test_run_pipeline_uses_the_default_gap_when_no_locus_declares_one(tmp_path, monkeypatch):
    assert _run_recording_gap(tmp_path, monkeypatch, ORDER_YML) == [25_000]


def test_run_pipeline_derives_the_gap_from_the_routed_locus(tmp_path, monkeypatch):
    """taxid=1 is in this order.yml's `taxonomic_scope`, so the route is
    `direct` and the locus's own curated gap is used."""
    assert _run_recording_gap(
        tmp_path, monkeypatch, WIDE_GAP_ORDER_YML, taxid=1
    ) == [50_000]


def test_run_pipeline_falls_back_to_the_default_gap_when_routing_failed(tmp_path, monkeypatch):
    """No taxid means `exhaustive`: the locus's wide gap is NOT inherited.

    Changed 2026-09-21 alongside the 120 kb Tremellales MAT gap. A run that
    could not route has no basis for clustering at the widest curated locus's
    distance -- see `family_registry.derive_max_cluster_gap`.
    """
    assert _run_recording_gap(tmp_path, monkeypatch, WIDE_GAP_ORDER_YML) == [
        25_000
    ]


def test_run_pipeline_takes_the_maximum_gap_over_routed_families(tmp_path, monkeypatch):
    """Two routed families, one default and one wide: the WIDER wins, because
    over-splitting a real locus is unrecoverable and under-splitting is not."""
    two_loci = WIDE_GAP_ORDER_YML + (
        "  - locus_name: bLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^b[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n"
    )
    assert _run_recording_gap(tmp_path, monkeypatch, two_loci, taxid=1) == [50_000]


def test_explicit_max_gap_argument_overrides_the_derived_value(tmp_path, monkeypatch):
    recorded = _run_recording_gap(tmp_path, monkeypatch, WIDE_GAP_ORDER_YML, max_gap=7_000)
    assert recorded == [7_000]


def test_evidence_diagnostics_record_rejected_clusters_not_just_admitted_ones(tmp_path):
    """END-TO-END guard on the diagnostics CALL SITE, not just the constant.

    `run_pipeline` must write a diagnostics row for EVERY (cluster, family)
    candidate -- every family with >=1 own hit in the cluster -- and mark the
    ones the real `evidence_floor` turned away with `admitted: false`. The
    rejected rows ARE the calibration dataset: they are the only record of what
    the floor is throwing away, and the next real rollout depends on them to
    replace this change's semi-synthetic measurement.

    This test fails if the diagnostics loop's `_DIAGNOSTICS_CANDIDATE_FLOOR` is
    ever "simplified" to the run's own `evidence_floor` (an easy "why are there
    two floors?" cleanup), because the rejected row would then never be
    written at all.

    Fixture: two clusters on two contigs. c1 has both of the family's genes and
    clears the default floor; c2 has one gene and is rejected by it.
    """
    import json

    _write_order(tmp_path)
    _write_record(tmp_path)
    diagnostics_path = tmp_path / "evidence_diagnostics.jsonl"

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None, **kwargs):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1",
                      "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                      "diamond_proteome"),
            # A lone gene on its own contig: 1 distinct gene, so the default
            # floor (>=2 genes) rejects this cluster.
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c2", 100, 200, "+", 31.0, "rec1",
                      "diamond_proteome"),
        ]

    run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa",
        taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=_no_localize,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        evidence_diagnostics_path=diagnostics_path,
    )

    rows = [json.loads(line) for line in diagnostics_path.read_text().splitlines()]
    by_contig = {row["contig"]: row for row in rows}
    assert set(by_contig) == {"c1", "c2"}, (
        f"both candidate clusters must be recorded, got {rows}"
    )
    assert by_contig["c1"]["admitted"] is True
    assert by_contig["c1"]["gene_count"] == 2
    # The load-bearing assertion: the REJECTED candidate is on disk.
    assert by_contig["c2"]["admitted"] is False
    assert by_contig["c2"]["gene_count"] == 1
    assert by_contig["c2"]["family"] == "P:aLocus"
