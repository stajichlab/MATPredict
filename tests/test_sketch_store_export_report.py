"""datasets/reads/sketch_store/export_report.py: the site report of a store; local samples counted but never named."""
import csv
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "datasets" / "reads" / "sketch_store"


def load(name):
    spec = importlib.util.spec_from_file_location(name, D / f"{name}.py"); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def write(path, cols, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def make_store(tmp):
    s = tmp / "store"; q = s / "query"
    ids = [f"SAMN{i}" for i in range(1, 9)] + ["local__secret_strain", "SAMN99"]
    src = {i: "local" if i.startswith("local__") else "public" for i in ids}
    calls = ["MAT1-1", "MAT1-2"] * 4 + ["MAT1-1", "MAT1-2"]
    write(s / "manifest.tsv", ["sample_id", "source", "run_accession", "strain", "study_accession", "first_public"],
          [dict(sample_id=i, source=src[i], run_accession="" if src[i] == "local" else "SRR" + i[4:], strain="st" + i[-1],
                study_accession="LabCollection" if src[i] == "local" else ("PRJNA9" if i == "SAMN99" else "PRJNA1")) for i in ids])
    (s / "sigs").mkdir()
    for i in ids:
        (s / "sigs" / f"{i}.sig.zip").write_bytes(b"x")
    write(q / "mat" / "mat_calls.tsv", ["sample_id", "source", "study_accession", "country", "sketch_call", "reads_type_call"],
          [dict(sample_id=i, source=src[i], study_accession="x", country="USA", sketch_call=c, reads_type_call=c) for i, c in zip(ids, calls)])
    write(q / "distances" / "samples_qc.tsv", ["sample_id", "mean_ani", "outlier", "clone_0.9995", "clone_0.998"],
          [dict(sample_id=i, mean_ani="0.83" if i == "SAMN99" else "0.995", outlier="1" if i == "SAMN99" else "0", **{"clone_0.9995": f"c{n}", "clone_0.998": "1"}) for n, i in enumerate(ids)])
    write(q / "distances" / "pcoa.tsv", ["sample_id", "PCo1", "PCo2"], [dict(sample_id=i, PCo1=str(n / 10), PCo2="0.1") for n, i in enumerate(ids)])
    write(q / "qc" / "qc.tsv", ["sample_id", "flag"], [])
    (q / "screen_t").mkdir()
    (q / "screen_t" / "summary.txt").write_text("bafA                                     len     638  specific hashes      3  present 4  partial 1  absent 5  | vs truth TP 2 FN 0 FP 1 TN 3\n")
    return s


def test_export_report(tmp_path):
    store = make_store(tmp_path)
    mf = tmp_path / "mf"
    load("mat_frequency").main(["--store", str(store), "--out", str(mf), "--min-group", "1"])
    out = tmp_path / "rep"
    load("export_report").main(["--store", str(store), "--matfreq", str(mf), "--species", "Species one", "--slug", "sp_one", "--out", str(out)])
    d = json.load(open(out / "population.json"))
    assert d["store"] == {**d["store"], "samples": 10, "public": 9, "local": 1, "sketched": 10, "missing": 0}
    assert d["excluded"] == {"kept": 9, "other_species": 1}
    assert d["overall"]["clone-corrected"]["n"] == 9
    assert d["outliers"]["kinds"] == {"other_species": 1} and d["outliers"]["other_species_studies"][0]["study"] == "PRJNA9"
    assert d["outliers"]["too_little_data_studies"] == []
    assert len(d["pcoa"]["points"]) == 9 and all(len(p) == 4 for p in d["pcoa"]["points"])       # the other species is not drawn
    assert d["screens"][0]["target"] == "bafA" and d["screens"][0]["truth"] == {"TP": 2, "FN": 0, "FP": 1, "TN": 3}
    text = (out / "population.json").read_text() + (out / "runs_public.tsv").read_text()
    assert "secret_strain" not in text                               # local isolates are counted, never named
    runs = list(csv.DictReader(open(out / "runs_public.tsv"), delimiter="\t"))
    assert len(runs) == 9 and {r["run"] for r in runs} >= {"SRR1", "SRR99"}
