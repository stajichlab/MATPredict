from __future__ import annotations
from pathlib import Path

from MATPredict.detect import pipeline as pipeline_module
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.polish import ExonSpan, PolishModel
from MATPredict.detect.search import SearchHit
from MATPredict.detect.pipeline import run_pipeline

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
    has a routed family with zero fast-path hits, so the zero-hit rescue cannot
    reach a real tblastn binary."""
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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


def test_fast_path_missing_gene_rescue_skips_localization_and_polishes_directly(tmp_path):
    """A fast-path cluster missing one core gene goes straight to
    polish_with_exonerate/polish_with_miniprot against a window padded around
    the EXISTING cluster's span, without ever calling search_localize."""
    _write_order(tmp_path)
    # A 100-aa curated mfa1 -> padding = 100 * 3 * 2.0 + 2000 = 2600 bp per side.
    _write_record(
        tmp_path,
        proteins=">rec1|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 100 + "\n",
    )
    localize_calls = []
    polish_calls = []

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                          "diamond_proteome")]

    def fake_localize(*args, **kwargs):
        localize_calls.append(args)
        return []

    def fake_exonerate(*, gene_name, window, **kwargs):
        polish_calls.append(("exonerate", gene_name, window))
        return _model("mfa1", "c1", 150, 260)

    def fake_miniprot(*, gene_name, window, **kwargs):
        polish_calls.append(("miniprot", gene_name, window))
        return _model("mfa1", "c1", 150, 260, method="miniprot_refine")

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_localize=fake_localize,
        polish_with_exonerate=fake_exonerate, polish_with_miniprot=fake_miniprot,
    )
    assert localize_calls == []  # Stage 1 skipped entirely
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


def test_polished_agree_and_disagree_produce_identical_tier(tmp_path):
    """Two otherwise-identical scenarios -- one where the two polish tools agree
    on a gene's boundaries and one where they disagree on the SAME gene -- must
    produce the same confidence tier and the same canonical coordinates.
    Agreement is reported, never consulted by assign_tier."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    coords = {"mfa1": (100, 200), "pra1": (300, 400)}

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return []

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        # aLocus has a foothold; bLocus has nothing at all.
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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


def test_rescued_cluster_with_an_unpolished_gene_is_capped_at_medium(tmp_path):
    """Polish eligibility is per-cluster, not gated on a global "genome-only"
    flag: a tblastn-rescued cluster on the FAST path is polished like any other
    localized cluster, so a gene neither tool can model there still caps the
    family's tier -- which a global flag would have silently prevented."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    coords = {"mfa1": (100, 200), "pra1": (300, 400)}

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return []

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_localize(genome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        # only pra1 is annotated; the cluster span is frozen at 300-400
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                          "diamond_proteome")]

    def polish(*, gene_name, **kwargs):
        return _model("mfa1", "c1", 150, 260)  # rescued gene, entirely left of 300

    outcome = run_pipeline(
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [SearchHit(three_gene_family.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                           "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )
    assert outcome.results[0].confidence == "low"


def test_genes_split_across_contigs_are_one_fragmented_multi_segment_call(tmp_path):
    """Finding 5: a family whose core genes land on different contigs, with no
    single cluster carrying them all, is one multi-segment locus with
    fragmented=True and a one-tier confidence downgrade."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\n" + "A" * 1000 + "\n>c2\n" + "A" * 1000 + "\n")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    outcome = run_pipeline(
        genome_fasta=genome, proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
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
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        ambiguity_floor=0.5,
    )
    assert len(outcome.results) == 2
    by_fragmented = {r.fragmented: r for r in outcome.results}
    assert True in by_fragmented and False in by_fragmented
    fragmented = by_fragmented[True]
    assert sorted({s.contig for s in fragmented.segments}) == ["c1", "c2"]
    independent = by_fragmented[False]
    assert independent.contig == "c3"
    assert sorted(independent.genes_found) == ["mfa1", "pra1"]


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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
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

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
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
