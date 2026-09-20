from __future__ import annotations
import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import (
    DetectionOutcome,
    DetectionResult,
    GeneEvidence,
    LocusSegment,
    NotDetectedFamily,
)
from MATPredict.detect.report import write_detection_gff3, write_detection_report

KEY = FamilyKey("Basidiomycota", "aLocus")

RESULT = DetectionResult(
    family_key=KEY, contig="c1", start=100, end=6443,
    confidence="high", idiomorph="undetermined", ambiguous_with=[],
    genes_found=["pra1"], genes_missing=["rba1"], fragmented=False,
    genes_not_searchable=["mfa1"],
    segments=[LocusSegment("c1", 100, 6443, contig_edge_distance=99)],
    gene_evidence=[
        GeneEvidence("pra1", "core_MAT", "c1", 1000, 2000, "+", 92.5, 87.0,
                     "5270_521_aLocus_a1", "diamond_proteome"),
    ],
    reference_records=["5270_521_aLocus_a1"],
)

OUTCOME = DetectionOutcome(
    results=[RESULT],
    not_detected=[
        NotDetectedFamily(
            family_key=FamilyKey("Basidiomycota", "bLocus"),
            reason="best cluster matched 0.25 of this family's expected genes, "
                   "below the ambiguity floor of 0.50",
            best_fraction_found=0.25,
            genes_found=["bE"],
            genes_missing=["bW"],
        )
    ],
    families_attempted=[KEY, FamilyKey("Basidiomycota", "bLocus")],
)


def test_write_detection_gff3_emits_locus_and_gene_features(tmp_path):
    out = tmp_path / "out.gff3"
    write_detection_gff3(OUTCOME, out)
    text = out.read_text()
    assert text.startswith("##gff-version 3")
    assert "##sequence-region c1 100 6443" in text
    assert "c1\tMATPredict\tMAT_locus\t100\t6443" in text
    # per-gene FEATURE lines, not just one locus-region line
    gene_line = next(line for line in text.splitlines() if "\tgene\t" in line and "Name=pra1" in line)
    assert "\t1000\t2000\t" in gene_line
    assert "role=core_MAT" in gene_line
    assert "present=true" in gene_line
    assert "identity=92.5" in gene_line
    assert "coverage=87.0" in gene_line
    assert "reference_record=5270_521_aLocus_a1" in gene_line
    # absent and not-searchable genes are explicit, with distinct semantics
    assert "Name=rba1;present=false" in text
    assert "Name=mfa1;present=false;not_searchable=true" in text


def test_write_detection_gff3_emits_one_sequence_region_per_segment(tmp_path):
    fragmented = DetectionResult(
        family_key=KEY, contig="c1", start=100, end=400,
        confidence="medium", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["mfa1", "pra1"], genes_missing=[], fragmented=True,
        segments=[LocusSegment("c1", 100, 200, 99), LocusSegment("c2", 300, 400, 299)],
    )
    out = tmp_path / "frag.gff3"
    write_detection_gff3(DetectionOutcome(results=[fragmented]), out)
    text = out.read_text()
    assert "##sequence-region c1 100 200" in text
    assert "##sequence-region c2 300 400" in text
    assert "fragmented=true" in text


def test_write_detection_gff3_gene_parent_is_scoped_to_its_own_contig(tmp_path):
    """Finding B (part 1) regression: a fragmented locus's gene on the second
    contig must NOT carry a Parent pointing at a MAT_locus feature declared
    only on the first contig -- that is not valid/clean GFF3 for a
    multi-contig feature set. Each segment gets its own MAT_locus feature,
    scoped to its own contig, and each gene's Parent points at the segment
    feature sharing its contig."""
    fragmented = DetectionResult(
        family_key=KEY, contig="c1", start=100, end=400,
        confidence="medium", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["mfa1", "pra1"], genes_missing=[], fragmented=True,
        segments=[LocusSegment("c1", 100, 200, 99), LocusSegment("c2", 300, 400, 299)],
        gene_evidence=[
            GeneEvidence("mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, None,
                         "5270_521_aLocus_a1", "diamond_proteome"),
            GeneEvidence("pra1", "core_MAT", "c2", 300, 400, "+", 95.0, None,
                         "5270_521_aLocus_a1", "diamond_proteome"),
        ],
    )
    out = tmp_path / "frag.gff3"
    write_detection_gff3(DetectionOutcome(results=[fragmented]), out)
    lines = out.read_text().splitlines()

    feature_lines = [line for line in lines if not line.startswith("#")]
    locus_lines = [line for line in feature_lines if "\tMAT_locus\t" in line]
    gene_lines = [line for line in feature_lines if "\tgene\t" in line]

    def _attr(line: str, key: str) -> str:
        attrs = dict(a.split("=", 1) for a in line.split("\t")[8].split(";"))
        return attrs[key]

    # one MAT_locus feature per segment, each declared on its own contig
    assert len(locus_lines) == 2
    locus_by_contig = {line.split("\t")[0]: line for line in locus_lines}
    c1_locus_id = _attr(locus_by_contig["c1"], "ID")
    c2_locus_id = _attr(locus_by_contig["c2"], "ID")
    assert c1_locus_id != c2_locus_id

    # every gene's Parent is on ITS OWN contig, never the other segment's contig
    mfa1_line = next(line for line in gene_lines if "Name=mfa1" in line)
    pra1_line = next(line for line in gene_lines if "Name=pra1" in line)
    assert mfa1_line.split("\t")[0] == "c1"
    assert _attr(mfa1_line, "Parent") == c1_locus_id
    assert pra1_line.split("\t")[0] == "c2"
    assert _attr(pra1_line, "Parent") == c2_locus_id


def test_write_detection_gff3_deduplicates_sequence_region_across_results(tmp_path):
    """Finding B (part 2) regression: two separate DetectionResults that both
    reference contig c1 must not each emit their own ##sequence-region c1
    pragma -- GFF3 tooling expects at most one per seqid."""
    key_b = FamilyKey("Basidiomycota", "bLocus")
    result_a = DetectionResult(
        family_key=KEY, contig="c1", start=100, end=200,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["mfa1"], genes_missing=[], fragmented=False,
    )
    result_b = DetectionResult(
        family_key=key_b, contig="c1", start=500, end=600,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["bE"], genes_missing=[], fragmented=False,
    )
    out = tmp_path / "dup.gff3"
    write_detection_gff3(DetectionOutcome(results=[result_a, result_b]), out)
    text = out.read_text()
    assert text.count("##sequence-region c1") == 1
    # the deduplicated pragma widens to cover both results' extents
    assert "##sequence-region c1 100 600" in text


def test_write_detection_report(tmp_path):
    out = tmp_path / "report.yaml"
    write_detection_report(OUTCOME, out)
    doc = yaml.safe_load(out.read_text())
    detected = doc["detected"][0]
    assert detected["family"] == "Basidiomycota:aLocus"
    assert detected["confidence"] == "high"
    assert detected["genes_missing"] == ["rba1"]
    assert detected["genes_not_searchable"] == ["mfa1"]
    assert detected["reference_records"] == ["5270_521_aLocus_a1"]
    assert detected["segments"][0]["contig_edge_distance"] == 99
    evidence = detected["gene_evidence"][0]
    assert evidence["gene"] == "pra1"
    assert evidence["identity"] == 92.5
    assert evidence["coverage"] == 87.0
    assert evidence["reference_record"] == "5270_521_aLocus_a1"


def test_write_detection_gff3_gene_feature_carries_status(tmp_path):
    out = tmp_path / "status.gff3"
    write_detection_gff3(OUTCOME, out)
    text = out.read_text()
    gene_line = next(line for line in text.splitlines() if "\tgene\t" in line and "Name=pra1" in line)
    assert "status=polished_agree" in gene_line


def test_write_detection_report_round_trips_disagree_status_and_alternate_model(tmp_path):
    """A polished_disagree gene's canonical coordinates AND its alternate
    model's data must both round-trip through the YAML report."""
    disagree_result = DetectionResult(
        family_key=KEY, contig="c1", start=100, end=6443,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["pra1"], genes_missing=[], fragmented=False,
        segments=[LocusSegment("c1", 100, 6443, contig_edge_distance=99)],
        gene_evidence=[
            GeneEvidence(
                "pra1", "core_MAT", "c1", 1000, 2000, "+", 92.5, 87.0,
                "5270_521_aLocus_a1", "exonerate_refine",
                status="polished_disagree",
                alternate_model={
                    "contig": "c1", "start": 1050, "end": 2050, "strand": "+",
                    "exons": [{"start": 1050, "end": 2050}],
                    "identity": 85.0, "method": "miniprot_refine",
                },
            ),
        ],
        reference_records=["5270_521_aLocus_a1"],
    )
    out = tmp_path / "disagree.yaml"
    write_detection_report(DetectionOutcome(results=[disagree_result]), out)
    doc = yaml.safe_load(out.read_text())
    evidence = doc["detected"][0]["gene_evidence"][0]
    # canonical data round-trips
    assert evidence["gene"] == "pra1"
    assert evidence["start"] == 1000
    assert evidence["end"] == 2000
    assert evidence["status"] == "polished_disagree"
    # alternate model's data round-trips too
    alt = evidence["alternate_model"]
    assert alt["contig"] == "c1"
    assert alt["start"] == 1050
    assert alt["end"] == 2050
    assert alt["method"] == "miniprot_refine"
    assert alt["identity"] == 85.0
    assert alt["exons"] == [{"start": 1050, "end": 2050}]


def test_write_detection_gff3_disagree_gene_carries_status_and_alt_attrs(tmp_path):
    disagree_result = DetectionResult(
        family_key=KEY, contig="c1", start=100, end=6443,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["pra1"], genes_missing=[], fragmented=False,
        segments=[LocusSegment("c1", 100, 6443, contig_edge_distance=99)],
        gene_evidence=[
            GeneEvidence(
                "pra1", "core_MAT", "c1", 1000, 2000, "+", 92.5, 87.0,
                "5270_521_aLocus_a1", "exonerate_refine",
                status="polished_disagree",
                alternate_model={
                    "contig": "c1", "start": 1050, "end": 2050, "strand": "+",
                    "exons": [{"start": 1050, "end": 2050}],
                    "identity": 85.0, "method": "miniprot_refine",
                },
            ),
        ],
    )
    out = tmp_path / "disagree.gff3"
    write_detection_gff3(DetectionOutcome(results=[disagree_result]), out)
    text = out.read_text()
    gene_line = next(line for line in text.splitlines() if "\tgene\t" in line and "Name=pra1" in line)
    assert "status=polished_disagree" in gene_line
    assert "alt_method=miniprot_refine" in gene_line
    assert "alt_start=1050" in gene_line
    assert "alt_end=2050" in gene_line


def test_not_polish_candidate_status_round_trips_through_yaml_and_gff3(tmp_path):
    """A gene found directly (never entered the polish pipeline) carries
    status=not_polish_candidate, distinct from status=unpolished, and that
    distinction must round-trip through both the YAML report and the GFF3
    gene-feature attributes."""
    result = DetectionResult(
        family_key=KEY, contig="c1", start=100, end=6443,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["pra1"], genes_missing=[], fragmented=False,
        segments=[LocusSegment("c1", 100, 6443, contig_edge_distance=99)],
        gene_evidence=[
            GeneEvidence(
                "pra1", "core_MAT", "c1", 1000, 2000, "+", 92.5, 87.0,
                "5270_521_aLocus_a1", "diamond_proteome",
                status="not_polish_candidate",
            ),
        ],
        reference_records=["5270_521_aLocus_a1"],
    )
    outcome = DetectionOutcome(results=[result])

    yaml_out = tmp_path / "not_polish_candidate.yaml"
    write_detection_report(outcome, yaml_out)
    doc = yaml.safe_load(yaml_out.read_text())
    evidence = doc["detected"][0]["gene_evidence"][0]
    assert evidence["status"] == "not_polish_candidate"
    assert evidence["status"] != "unpolished"
    assert evidence["alternate_model"] is None

    gff3_out = tmp_path / "not_polish_candidate.gff3"
    write_detection_gff3(outcome, gff3_out)
    text = gff3_out.read_text()
    gene_line = next(line for line in text.splitlines() if "\tgene\t" in line and "Name=pra1" in line)
    assert "status=not_polish_candidate" in gene_line
    assert "status=unpolished" not in gene_line


def _write_genome_fasta(tmp_path, contig_id, prefix_len, orf, suffix_len):
    """A real synthetic single-contig genome FASTA with a real ORF starting
    right after `prefix_len` bases (so the ORF's 1-based start is
    `prefix_len + 1`), for `_extract_translated_gene` to genuinely slice
    and translate."""
    sequence = ("A" * prefix_len) + orf + ("A" * suffix_len)
    path = tmp_path / "genome.fasta"
    path.write_text(f">{contig_id}\n{sequence}\n")
    return path


# ATG + 10x GCT (Ala) + TAA stop -- translates to "MAAAAAAAAAA" (stop dropped).
_ORF = "ATG" + "GCT" * 10 + "TAA"
_ORF_START = 100  # 1-based
_ORF_END = _ORF_START + len(_ORF) - 1  # 135
_ORF_PROTEIN = "MAAAAAAAAAA"


def test_write_detection_gff3_with_genome_fasta_writes_companion_fasta_and_cds(tmp_path):
    genome_fasta = _write_genome_fasta(
        tmp_path, "c1", prefix_len=_ORF_START - 1, orf=_ORF, suffix_len=20,
    )
    result = DetectionResult(
        family_key=KEY, contig="c1", start=_ORF_START, end=_ORF_END,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["pra1"], genes_missing=[], fragmented=False,
        segments=[LocusSegment("c1", _ORF_START, _ORF_END, contig_edge_distance=99)],
        gene_evidence=[
            GeneEvidence("pra1", "core_MAT", "c1", _ORF_START, _ORF_END, "+", 92.5, 87.0,
                         "5270_521_aLocus_a1", "diamond_proteome"),
        ],
        reference_records=["5270_521_aLocus_a1"],
    )
    outcome = DetectionOutcome(results=[result], not_detected=[], families_attempted=[KEY])

    out = tmp_path / "with_seq.gff3"
    write_detection_gff3(outcome, out, genome_fasta=genome_fasta)

    # (a) companion FASTA written with the real sliced sequence, same contig name.
    fasta_out = out.with_suffix(".fasta")
    assert fasta_out.exists()
    fasta_text = fasta_out.read_text()
    assert fasta_text.startswith(">c1\n")
    # The sequence is wrapped at 60 columns, so join the sequence lines back
    # up before looking for the ORF.
    sequence = "".join(
        line for line in fasta_text.splitlines() if not line.startswith(">")
    )
    assert _ORF in sequence

    # (b) a CDS feature line with a translation= attribute for the gene.
    text = out.read_text()
    cds_line = next(line for line in text.splitlines() if "\tCDS\t" in line)
    assert f"\t{_ORF_START}\t{_ORF_END}\t" in cds_line
    assert f"translation={_ORF_PROTEIN}" in cds_line
    assert "Parent=" in cds_line
    # the gene feature is still emitted, unaffected by the new CDS feature.
    gene_line = next(line for line in text.splitlines() if "\tgene\t" in line and "Name=pra1" in line)
    assert gene_line != cds_line


def test_write_detection_gff3_contig_mismatch_falls_back_gracefully(tmp_path):
    """A genome FASTA that doesn't contain the report's contig must not crash
    the write -- no companion FASTA, no CDS feature, gene feature still
    written normally."""
    genome_fasta = _write_genome_fasta(
        tmp_path, "other_contig", prefix_len=_ORF_START - 1, orf=_ORF, suffix_len=20,
    )
    out = tmp_path / "mismatch.gff3"
    write_detection_gff3(OUTCOME, out, genome_fasta=genome_fasta)

    text = out.read_text()
    assert "\tCDS\t" not in text
    gene_line = next(line for line in text.splitlines() if "\tgene\t" in line and "Name=pra1" in line)
    assert gene_line  # gene feature still present
    assert not out.with_suffix(".fasta").exists()


def test_write_detection_gff3_without_genome_fasta_is_unchanged(tmp_path):
    """Omitting `genome_fasta` (the default) must produce byte-identical
    output to before this parameter existed -- no regression for existing
    callers."""
    baseline = tmp_path / "baseline.gff3"
    write_detection_gff3(OUTCOME, baseline)

    explicit_none = tmp_path / "explicit_none.gff3"
    write_detection_gff3(OUTCOME, explicit_none, genome_fasta=None)

    assert baseline.read_bytes() == explicit_none.read_bytes()
    assert not baseline.with_suffix(".fasta").exists()
    assert not explicit_none.with_suffix(".fasta").exists()


def test_write_detection_report_lists_not_detected_families(tmp_path):
    """Sub-floor families must appear with a reason, never be silently dropped."""
    out = tmp_path / "report.yaml"
    write_detection_report(OUTCOME, out)
    doc = yaml.safe_load(out.read_text())
    assert doc["families_attempted"] == ["Basidiomycota:aLocus", "Basidiomycota:bLocus"]
    assert doc["not_detected"][0]["family"] == "Basidiomycota:bLocus"
    assert "below the ambiguity floor" in doc["not_detected"][0]["reason"]
    assert doc["not_detected"][0]["best_fraction_found"] == 0.25


def test_write_detection_gff3_parses_the_genome_exactly_once(tmp_path, monkeypatch):
    """The genome FASTA must be opened/parsed ONCE per call, not once per
    gene plus once more for the companion FASTA.

    Before this test, `write_detection_gff3` called
    `_extract_translated_gene(genome_fasta, ...)` inside the per-gene loop --
    and that helper opens and `SeqIO.parse`s the WHOLE genome on every call --
    then opened the genome a final time to build the companion FASTA. For a
    real fungal genome (tens of MB) and a locus with several genes that is a
    full re-parse per gene, entirely avoidable because every gene's sequence
    comes from the same small set of contigs the companion FASTA already
    reads. The assertion counts real calls to the module's own FASTA-opening
    helper rather than measuring elapsed time, so it cannot pass or fail for
    reasons unrelated to the number of parses.
    """
    from Bio import SeqIO

    genome_fasta = _write_genome_fasta(
        tmp_path, "c1", prefix_len=_ORF_START - 1, orf=_ORF, suffix_len=20,
    )
    # `SeqIO.parse` is counted rather than either module's own
    # `_open_fasta_text`, because a parse could be issued from `report.py` or
    # from `benchmark.py`'s `_extract_translated_gene` (which binds its own
    # module-level helper); counting the single function both paths must go
    # through catches the N+1 wherever it lives.
    real_parse = SeqIO.parse
    parses = []

    def counting_parse(handle, fmt, *args, **kwargs):
        parses.append(fmt)
        return real_parse(handle, fmt, *args, **kwargs)

    monkeypatch.setattr(SeqIO, "parse", counting_parse)

    # Two genes on the same contig -- the N+1 shape this test exists to catch
    # only shows up with more than one gene.
    evidence = [
        GeneEvidence("pra1", "core_MAT", "c1", _ORF_START, _ORF_END, "+", 92.5, 87.0,
                     "5270_521_aLocus_a1", "diamond_proteome"),
        GeneEvidence("rba1", "core_MAT", "c1", _ORF_START, _ORF_END, "+", 90.0, 85.0,
                     "5270_521_aLocus_a1", "diamond_proteome"),
    ]
    result = DetectionResult(
        family_key=KEY, contig="c1", start=_ORF_START, end=_ORF_END,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["pra1", "rba1"], genes_missing=[], fragmented=False,
        segments=[LocusSegment("c1", _ORF_START, _ORF_END, contig_edge_distance=99)],
        gene_evidence=evidence,
        reference_records=["5270_521_aLocus_a1"],
    )
    outcome = DetectionOutcome(results=[result], not_detected=[], families_attempted=[KEY])

    out = tmp_path / "once.gff3"
    write_detection_gff3(outcome, out, genome_fasta=genome_fasta)

    assert len(parses) == 1, f"genome parsed {len(parses)} times, expected exactly 1"
    # ...and the output is still correct: both genes got a real CDS feature.
    cds_lines = [line for line in out.read_text().splitlines() if "\tCDS\t" in line]
    assert len(cds_lines) == 2
    assert all(f"translation={_ORF_PROTEIN}" in line for line in cds_lines)


def test_companion_fasta_is_wrapped_at_60_columns(tmp_path):
    """Sequence lines must be wrapped at the standard 60-column FASTA width.

    A real fungal contig is megabases long; writing it as one line makes the
    companion FASTA unreadable in a pager and is rejected or mangled by some
    downstream parsers. 60 is the width the rest of this project emits.
    """
    # A contig comfortably longer than 60 bases, so wrapping is observable.
    genome_fasta = _write_genome_fasta(
        tmp_path, "c1", prefix_len=_ORF_START - 1, orf=_ORF, suffix_len=200,
    )
    result = DetectionResult(
        family_key=KEY, contig="c1", start=_ORF_START, end=_ORF_END,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["pra1"], genes_missing=[], fragmented=False,
        segments=[LocusSegment("c1", _ORF_START, _ORF_END, contig_edge_distance=99)],
        gene_evidence=[
            GeneEvidence("pra1", "core_MAT", "c1", _ORF_START, _ORF_END, "+", 92.5, 87.0,
                         "5270_521_aLocus_a1", "diamond_proteome"),
        ],
        reference_records=["5270_521_aLocus_a1"],
    )
    outcome = DetectionOutcome(results=[result], not_detected=[], families_attempted=[KEY])

    out = tmp_path / "wrapped.gff3"
    write_detection_gff3(outcome, out, genome_fasta=genome_fasta)

    lines = out.with_suffix(".fasta").read_text().splitlines()
    assert lines[0] == ">c1"
    sequence_lines = [line for line in lines if not line.startswith(">")]
    assert sequence_lines, "companion FASTA has no sequence lines"
    assert all(len(line) <= 60 for line in sequence_lines)
    assert len(sequence_lines) > 1, "sequence was not wrapped at all"
    # Wrapping must not alter the sequence itself.
    assert _ORF in "".join(sequence_lines)


def test_write_detection_report_records_the_routing_decision(tmp_path):
    """Task 1 item 4 / the plan's global constraint: every routing fallback
    must be visible in the report, never silent. A reader must be able to
    tell "these 2 families were searched because the taxid's phylum was
    Ascomycota" apart from "these 2 were searched because their scope
    actually matched"."""
    out = tmp_path / "routed.yaml"
    outcome = DetectionOutcome(
        results=[],
        not_detected=[],
        families_attempted=[FamilyKey("Ascomycota", "MATsc"), FamilyKey("Ascomycota", "MATyl")],
        routing_mode="phylum_fallback",
    )
    write_detection_report(outcome, out)
    doc = yaml.safe_load(out.read_text())
    assert doc["routing_mode"] == "phylum_fallback"
    assert doc["families_attempted"] == ["Ascomycota:MATsc", "Ascomycota:MATyl"]


def test_write_detection_report_routing_mode_defaults_to_null(tmp_path):
    """An outcome built without routing information (a direct `run_pipeline`
    call in a test, say) still writes the key, as an explicit null, rather
    than omitting it and making a consumer guess."""
    out = tmp_path / "unrouted.yaml"
    write_detection_report(OUTCOME, out)
    assert yaml.safe_load(out.read_text())["routing_mode"] is None
