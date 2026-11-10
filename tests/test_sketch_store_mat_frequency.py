"""datasets/reads/sketch_store/mat_frequency.py: call choice, exclusions, clone correction, Wilson interval and binomial test."""
import csv
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("mat_frequency", ROOT / "datasets" / "reads" / "sketch_store" / "mat_frequency.py")
mf = importlib.util.module_from_spec(spec); spec.loader.exec_module(mf)


def write(path, cols, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def store(tmp_path):
    s = tmp_path / "store"; q = s / "query"
    # c1: a clone sequenced 4 times (MAT1-1); s5 MAT1-2 alone; s6 reads-type low_depth but sketch MAT1-2; s7 other species;
    # s8 excluded by the curator; s9 QC mixture; s10 'both'
    mat = [dict(sample_id=f"s{i}", source="public", study_accession="P1" if i < 5 else "P2", country="USA: CA" if i < 7 else "",
                sketch_call="MAT1-1", reads_type_call="MAT1-1") for i in range(1, 5)]
    mat += [dict(sample_id="s5", source="public", study_accession="P2", sketch_call="MAT1-2", reads_type_call="MAT1-2"),
            dict(sample_id="s6", source="local", study_accession="L", sketch_call="MAT1-2", reads_type_call="low_depth"),
            dict(sample_id="s7", source="public", study_accession="P3", sketch_call="MAT1-1", reads_type_call="MAT1-1"),
            dict(sample_id="s8", source="local", study_accession="L", sketch_call="MAT1-1", reads_type_call="MAT1-1"),
            dict(sample_id="s9", source="public", study_accession="P2", sketch_call="MAT1-2", reads_type_call="MAT1-2"),
            dict(sample_id="s10", source="public", study_accession="P2", sketch_call="both", reads_type_call="both")]
    write(q / "mat" / "mat_calls.tsv", ["sample_id", "source", "study_accession", "country", "sketch_call", "reads_type_call"], mat)
    sq = [dict(sample_id=f"s{i}", mean_ani="0.995", **{"clone_0.9995": "c1", "clone_0.998": "L1"}) for i in range(1, 5)]
    sq += [dict(sample_id=s, mean_ani="0.995", **{"clone_0.9995": c, "clone_0.998": "L2"}) for s, c in
           (("s5", "c5"), ("s6", "c6"), ("s8", "c8"), ("s9", "c9"), ("s10", "c10"))]
    sq += [dict(sample_id="s7", mean_ani="0.83")]
    write(q / "distances" / "samples_qc.tsv", ["sample_id", "mean_ani", "clone_0.9995", "clone_0.998"], sq)
    write(q / "qc" / "qc.tsv", ["sample_id", "flag"], [dict(sample_id="s9", flag="mixture?")])
    write(s / "exclude.tsv", ["sample_id", "reason"], [dict(sample_id="s8", reason="low quality")])
    write(s / "manifest.tsv", ["sample_id", "strain"], [dict(sample_id=f"s{i}", strain=f"st{i}") for i in range(1, 11)])
    return s


def test_mat_frequency(tmp_path):
    out = tmp_path / "out"
    mf.main(["--store", str(store(tmp_path)), "--out", str(out), "--min-group", "1"])
    samples = {r["sample_id"]: r for r in csv.DictReader(open(out / "samples.tsv"), delimiter="\t")}
    assert samples["s6"]["call"] == "MAT1-2"                      # reads-type low depth -> sketch call
    assert samples["s7"]["excluded"] == "other_species"
    assert samples["s8"]["excluded"] == "curator_exclude"
    assert samples["s9"]["excluded"] == "qc_mixture"
    assert samples["s10"]["call"] == "both"
    assert samples["s1"]["country"] == "USA"
    per = next(csv.DictReader(open(out / "freq_overall_samples.tsv"), delimiter="\t"))
    assert (per["MAT1_1"], per["MAT1_2"], per["both"]) == ("4", "2", "1")
    cl = next(csv.DictReader(open(out / "freq_overall_clone-corrected.tsv"), delimiter="\t"))
    assert (cl["MAT1_1"], cl["MAT1_2"]) == ("1", "2")              # the 4-sample clone votes once
    assert (out / "freq_lineage_clones.tsv").exists() and "clone groups" in (out / "summary.txt").read_text()


def test_stats():
    lo, hi = mf.wilson(58, 100)
    assert 0.48 < lo < 0.49 and 0.67 < hi < 0.68
    assert abs(mf.binom_p(5, 10) - 1.0) < 1e-9
    assert abs(mf.binom_p(0, 10) - 2 / 1024) < 1e-12
