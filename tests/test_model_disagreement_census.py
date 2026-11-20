"""Fixture test of scripts/model_disagreement_census.py (#290, #349): classes, trigger, span change and the tar stream."""
from __future__ import annotations

import gzip
import importlib.util
import io
import sys
import tarfile
from pathlib import Path

import yaml

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("mdcen", SCRIPTS / "model_disagreement_census.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

GENOME = "CCCC" + "ATGGCTGCTAAA" + "GT" + "A" * 40 + "AG" + "GCTAAAGGTTAA" + "CCCC"
GOOD = [{"start": 5, "end": 16}, {"start": 61, "end": 72}]
SHIFTED = [{"start": 5, "end": 18}, {"start": 63, "end": 72}]


def entry(primary, alternate):
    return {"gene": "G1", "role": "core_MAT", "status": "polished_disagree", "contig": "c1", "start": 5, "end": 72, "strand": "+", "method": "exonerate_refine",
            "identity": 99.0, "exons": primary, "alternate_model": {"contig": "c1", "start": 5, "end": 72, "strand": "+", "method": "miniprot_refine", "identity": 90.0, "exons": alternate}}


def doc(*entries, code=1):
    return {"genetic_code": code, "detected": [{"family": "Ascomycota:MAT", "confidence": "high", "start": 5, "end": 72, "gene_evidence": list(entries)}]}


def test_classify():
    assert mod.classify([], []) == "both_clean" and mod.classify(["x"], []) == "primary_issue_alt_clean"
    assert mod.classify([], ["x"]) == "primary_clean_alt_issue" and mod.classify(["x"], ["y"]) == "both_issue"


def test_run_classifies_and_flags_the_noncanonical_trigger(tmp_path):
    (tmp_path / "g").mkdir()
    for s in ("A", "B", "C"):
        with gzip.open(tmp_path / "g" / f"{s}_x.fasta.gz", "wt") as fh:
            fh.write(">c1\n" + GENOME + "\n")
    docs = [("A", doc(entry(SHIFTED, GOOD))), ("B", doc(entry(GOOD, SHIFTED))), ("C", doc(entry(GOOD, GOOD))), ("D", doc(entry(GOOD, GOOD)))]
    rows = mod.run(docs, tmp_path / "g", tmp_path / "db", tmp_path / "out" / "t", "t")
    by = {r["sample"]: r for r in rows}
    assert by["A"]["class"] == "primary_issue_alt_clean" and by["A"]["trigger_noncanonical"] is True and by["A"]["d_aa"] == 0 and by["A"]["d_end"] == 0
    assert by["B"]["class"] == "primary_clean_alt_issue" and not by["B"]["trigger_noncanonical"]
    assert by["C"]["class"] == "both_clean"
    assert by["D"]["class"] == "not_checked"                 # no genome file for D
    summ = (tmp_path / "out" / "t.summary.tsv").read_text()
    assert "noncanonical_primary_alt_none\tgenes\t1" in summ and "span_changes" in summ


def test_unsupported_genetic_code_is_not_checked(tmp_path):
    (tmp_path / "g").mkdir()
    with gzip.open(tmp_path / "g" / "A_x.fasta.gz", "wt") as fh:
        fh.write(">c1\n" + GENOME + "\n")
    rows = mod.run([("A", doc(entry(GOOD, GOOD), code=4))], tmp_path / "g", tmp_path / "db", tmp_path / "o", "t")
    assert rows[0]["class"] == "not_checked"


def test_tar_stream_reader():
    buf = io.BytesIO()
    payload = yaml.safe_dump(doc(entry(GOOD, GOOD))).encode()
    with tarfile.open(fileobj=buf, mode="w") as tf:
        info = tarfile.TarInfo("wave_0/small_000/runs/S1/detection_report.yaml")
        info.size = len(payload)
        tf.addfile(info, io.BytesIO(payload))
        other = tarfile.TarInfo("wave_0/small_000/runs/S1/wall_seconds")
        other.size = 2
        tf.addfile(other, io.BytesIO(b"1\n"))
    buf.seek(0)
    got = list(mod._iter_tar(buf))
    assert [g[0] for g in got] == ["S1"] and got[0][1]["genetic_code"] == 1


def test_inputs_tsv_maps_prefixed_samples_to_genome_paths(tmp_path):
    fa = tmp_path / "somewhere" / "G.sorted.fasta"
    fa.parent.mkdir()
    fa.write_text(">c1\n" + GENOME + "\n")
    (tmp_path / "inputs.tsv").write_text(f"source\tid\tgenome\torder\nLCG\tG1\t{fa}\tMucorales\n")
    inputs = mod.read_inputs(tmp_path / "inputs.tsv")
    assert inputs == {"LCG__G1": fa}
    rows = mod.run([("LCG__G1", doc(entry(GOOD, GOOD)))], tmp_path / "none", tmp_path / "db", tmp_path / "o", "t", inputs)
    assert rows[0]["class"] == "both_clean"
