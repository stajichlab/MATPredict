"""The pre-sign-off regression check (curator ruling 2026-09-30).

Before a new record or classifier rebuild is signed off, every existing call
whose presence, label, confidence, locus class, classifier verdict, gene set or
core gene model changed must be listed -- no silent aggregation.
"""
import json

import yaml

from MATPredict.detect.regression import compare_runs, load_reports, write_outputs


def _call(contig="c1", start=100, end=5000, idiomorph="Plus", confidence="high",
          locus_class="mat_locus", clf_input="model", scores=(30.0, 110.0), margin=80.0,
          genes=("sexP", "rnhA"), exons=((1000, 1900),), identity=60.0):
    return {
        "family": "Mucoromycota:MAT", "contig": contig, "start": start, "end": end,
        "idiomorph": idiomorph, "confidence": confidence, "locus_class": locus_class,
        "verification": None, "genes_found": list(genes),
        "idiomorph_classifier": {"classifier_input": clf_input, "margin": margin,
                                 "scores": {"Minus": scores[0], "Plus": scores[1]}},
        "gene_evidence": [
            {"gene": "sexP", "role": "core_MAT", "contig": contig, "start": exons[0][0],
             "end": exons[-1][1], "identity": identity,
             "exons": [{"start": s, "end": e} for s, e in exons]},
            {"gene": "rnhA", "role": "flanking_conserved", "contig": contig, "start": 3000,
             "end": 4000, "identity": 90.0, "exons": [{"start": 3000, "end": 4000}]},
        ],
    }


def _write(root, genome, detected=(), suppressed=(), diagnostics=()):
    d = root / "runs" / genome
    d.mkdir(parents=True)
    (d / "detection_report.yaml").write_text(yaml.safe_dump(
        {"detected": list(detected), "suppressed_loci": list(suppressed)}))
    if diagnostics:
        (d / "evidence_diagnostics.jsonl").write_text(
            "\n".join(json.dumps(x) for x in diagnostics) + "\n")


def _changes(rows, genome):
    return {c for r in rows if r["genome"] == genome for c in r["change_types"].split(",") if c}


def test_identical_runs_report_no_change(tmp_path):
    for side in ("base", "cand"):
        _write(tmp_path / side, "g1", detected=[_call()])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    assert len(rows) == 1 and rows[0]["change_types"] == ""


def test_a_lost_call_names_the_candidate_withheld_reason(tmp_path):
    _write(tmp_path / "base", "g1", detected=[_call()])
    _write(tmp_path / "cand", "g1", suppressed=[{
        "family": "Mucoromycota:MAT", "contig": "c1", "start": 200, "end": 4800,
        "idiomorph": "undetermined", "genes_found": ["sexP", "rnhA"],
        "withheld_reason": "mat_gene_gate", "idiomorph_classifier": None}])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    assert "call_lost" in _changes(rows, "g1")
    assert rows[0]["cand_state"] == "withheld:mat_gene_gate"


def test_a_gained_call_is_listed(tmp_path):
    _write(tmp_path / "base", "g1")
    _write(tmp_path / "cand", "g1", detected=[_call(idiomorph="Minus")])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    assert _changes(rows, "g1") == {"call_gained"}
    assert rows[0]["base_state"] == "absent"


def test_label_confidence_and_class_changes_are_each_flagged(tmp_path):
    _write(tmp_path / "base", "g1", detected=[_call()])
    _write(tmp_path / "cand", "g1", detected=[_call(
        idiomorph="undetermined", confidence="medium", locus_class="partial_locus")])
    ch = _changes(compare_runs(load_reports(tmp_path / "base"),
                               load_reports(tmp_path / "cand")), "g1")
    assert {"idiomorph_changed", "confidence_changed", "locus_class_changed"} <= ch


def test_a_changed_core_model_is_flagged_even_when_the_label_holds(tmp_path):
    # The Circinella case: same label, but the sexP model detection built changed.
    _write(tmp_path / "base", "g1", detected=[_call(exons=((1000, 1900),), identity=60.0)])
    _write(tmp_path / "cand", "g1", detected=[_call(exons=((1000, 1600),), identity=52.0)])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    assert "core_model_changed" in _changes(rows, "g1")
    assert "sexP" in rows[0]["model_changes"]


def test_classifier_shift_below_threshold_is_not_a_change_but_is_recorded(tmp_path):
    _write(tmp_path / "base", "g1", detected=[_call(margin=80.0)])
    _write(tmp_path / "cand", "g1", detected=[_call(margin=82.0)])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"),
                        score_delta=5.0)
    assert rows[0]["change_types"] == ""
    assert rows[0]["margin_delta"] == 2.0


def test_classifier_shift_above_threshold_and_input_change_are_flagged(tmp_path):
    _write(tmp_path / "base", "g1", detected=[_call(margin=80.0, clf_input="model")])
    _write(tmp_path / "cand", "g1", detected=[_call(margin=60.0, clf_input="hsp_fragment")])
    ch = _changes(compare_runs(load_reports(tmp_path / "base"),
                               load_reports(tmp_path / "cand"), score_delta=5.0), "g1")
    assert {"classifier_shift", "classifier_input_changed"} <= ch


def test_gene_set_change_is_flagged(tmp_path):
    _write(tmp_path / "base", "g1", detected=[_call(genes=("sexP", "rnhA"))])
    _write(tmp_path / "cand", "g1", detected=[_call(genes=("sexP", "rnhA", "glrA"))])
    ch = _changes(compare_runs(load_reports(tmp_path / "base"),
                               load_reports(tmp_path / "cand")), "g1")
    assert "gene_set_changed" in ch


def test_withheld_reason_marks_polish_cap_from_diagnostics(tmp_path):
    # The M. pusillus case: the cluster was skipped by the polish cap.
    _write(tmp_path / "base", "g1", detected=[_call()])
    _write(tmp_path / "cand", "g1",
           suppressed=[{"family": "Mucoromycota:MAT", "contig": "c1", "start": 100,
                        "end": 5000, "idiomorph": "undetermined", "genes_found": ["sexP"],
                        "withheld_reason": "modelled_gene_bar"}],
           diagnostics=[{"kind": "evidence", "family": "Mucoromycota:MAT", "contig": "c1",
                         "cluster_start": 100, "cluster_end": 5000, "polish_capped": True}])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    assert rows[0]["cand_state"] == "withheld:modelled_gene_bar+polish_capped"


def test_polish_cap_read_from_zstd_compressed_diagnostics(tmp_path):
    # Archived runs keep evidence_diagnostics.jsonl compressed as .jsonl.zst.
    from MATPredict.db.local_taxonomy import _zstd_compress
    _write(tmp_path / "cand", "g1",
           diagnostics=[{"kind": "evidence", "family": "Mucoromycota:MAT", "contig": "c1",
                         "cluster_start": 100, "cluster_end": 5000, "polish_capped": True}])
    plain = tmp_path / "cand" / "runs" / "g1" / "evidence_diagnostics.jsonl"
    plain.with_suffix(".jsonl.zst").write_bytes(
        _zstd_compress(plain.read_bytes()))
    plain.unlink()
    rep = load_reports(tmp_path / "cand")["g1"]
    assert rep["_capped"] == [("Mucoromycota:MAT", "c1", 100, 5000)]


def test_loci_match_by_contig_overlap_not_exact_coordinates(tmp_path):
    _write(tmp_path / "base", "g1", detected=[_call(start=100, end=5000)])
    _write(tmp_path / "cand", "g1", detected=[_call(start=90, end=5600)])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    assert len(rows) == 1 and "span_changed" in _changes(rows, "g1")


def test_genome_present_on_one_side_only_is_reported(tmp_path):
    _write(tmp_path / "base", "g1", detected=[_call()])
    (tmp_path / "cand").mkdir()
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    assert _changes(rows, "g1") == {"genome_missing_in_candidate"}


def test_nested_input_directories_keep_distinct_keys(tmp_path):
    # Zygo runs: <out>/scaffold/runs/<org> and <out>/contig/runs/<org>.
    for inp in ("scaffold", "contig"):
        d = tmp_path / inp / "runs" / "orgA"
        d.mkdir(parents=True)
        (d / "detection_report.yaml").write_text(yaml.safe_dump({"detected": [_call()]}))
    assert set(load_reports(tmp_path)) == {"scaffold/orgA", "contig/orgA"}


def test_outputs_list_every_changed_locus(tmp_path):
    _write(tmp_path / "base", "g1", detected=[_call()])
    _write(tmp_path / "base", "g2", detected=[_call()])
    _write(tmp_path / "cand", "g1", detected=[_call(idiomorph="Minus")])
    _write(tmp_path / "cand", "g2", detected=[_call()])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    tsv, md = write_outputs(rows, tmp_path / "out", title="t", score_delta=5.0)
    text = md.read_text()
    assert "g1" in text and "idiomorph_changed" in text
    assert "g2" not in text.split("## Changes to calls")[1]
    assert len(tsv.read_text().strip().splitlines()) == 3


def test_withheld_only_changes_go_to_the_appendix_not_the_main_list(tmp_path):
    w = {"family": "Mucoromycota:MAT", "contig": "c9", "start": 1, "end": 900,
         "idiomorph": "undetermined", "genes_found": ["tptA"], "withheld_reason": "modelled_gene_bar"}
    _write(tmp_path / "base", "g1", detected=[_call()], suppressed=[w])
    _write(tmp_path / "cand", "g1", detected=[_call(idiomorph="Minus")],
           suppressed=[{**w, "withheld_reason": "below_fraction_floor"}])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    _, md = write_outputs(rows, tmp_path / "out", title="t")
    main = md.read_text()
    appendix = (tmp_path / "out" / "regression_withheld_changes.md").read_text()
    assert "c1" in main and "c9" not in main.split("## Changes to calls")[1]
    assert "c9" in appendix and "withheld_reason_changed" in appendix


def test_a_lost_call_does_not_claim_its_genes_were_removed(tmp_path):
    _write(tmp_path / "base", "g1", detected=[_call()])
    _write(tmp_path / "cand", "g1", suppressed=[{
        "family": "Mucoromycota:MAT", "contig": "c1", "start": 100, "end": 5000,
        "idiomorph": "Plus", "genes_found": ["sexP", "rnhA"], "withheld_reason": "mat_gene_gate",
        "idiomorph_classifier": {"classifier_input": "model", "margin": 54.2,
                                 "scores": {"Minus": 32.2, "Plus": 86.4}}}])
    rows = compare_runs(load_reports(tmp_path / "base"), load_reports(tmp_path / "cand"))
    r = rows[0]
    assert r["model_changes"] == "no_gene_evidence_in_candidate"
    assert "core_model_changed" not in r["change_types"]
    assert (r["base_best_score"], r["cand_best_score"]) == (110.0, 86.4)


def test_cli_diffs_named_pairs_and_writes_an_index(tmp_path):
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location(
        "regression_check", Path(__file__).resolve().parents[2] / "scripts" / "regression_check.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for side in ("b", "c"):
        _write(tmp_path / side, "g1", detected=[_call(idiomorph="Plus" if side == "b" else "Minus")])
    rc = mod.main(["diff", "--out", str(tmp_path / "out"), "--title", "t",
                   "--pair", "p1", str(tmp_path / "b"), str(tmp_path / "c"),
                   "--pair", "p2", str(tmp_path / "b"), str(tmp_path / "b")])
    assert rc == 0
    index = (tmp_path / "out" / "summary.md").read_text()
    assert "| p1 | 1 | 1 |" in index and "| p2 | 1 | 0 |" in index
