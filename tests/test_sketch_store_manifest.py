"""datasets/reads/sketch_store/build_manifest.py: one run per BioSample, flagged studies dropped, 30x target spots, local rows."""
import csv
import gzip
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "datasets" / "reads" / "sketch_store" / "build_manifest.py"
COLS = ["run_accession", "study_accession", "sample_accession", "scientific_name", "tax_id", "instrument_platform", "instrument_model",
        "library_layout", "library_strategy", "library_source", "read_count", "base_count", "first_public", "center_name", "country",
        "collection_date", "strain", "isolate", "sample_title", "study_title", "fastq_ftp", "fastq_aspera"]


def run(tmp_path, runs, local=None):
    rp = tmp_path / "runs.tsv.gz"
    with gzip.open(rp, "wt", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t"); w.writeheader()
        for r in runs:
            w.writerow({c: r.get(c, "") for c in COLS})
    tp = tmp_path / "tax.csv.gz"
    with gzip.open(tp, "wt") as fh:
        fh.write("taxid,species_taxid\n746128,746128\n330879,746128\n5062,5062\n")
    out = tmp_path / "manifest.tsv"
    cmd = [sys.executable, str(SCRIPT), "--runs", str(rp), "--taxmap", str(tp), "--species-taxid", "746128", "--genome-mb", "30", "--out", str(out)]
    if local:
        lp = tmp_path / "samples.csv"; lp.write_text(local); cmd += ["--local-samples", str(lp)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return list(csv.DictReader(open(out), delimiter="\t"))


def r(run, bs, study, layout="PAIRED", reads=10_000_000, bases=3_000_000_000, tax="746128", strain="", title="Global isolates"):
    return dict(run_accession=run, sample_accession=bs, study_accession=study, library_layout=layout, read_count=reads, base_count=bases,
                tax_id=tax, strain=strain, study_title=title, scientific_name="Aspergillus fumigatus")


def test_one_run_per_biosample_paired_and_largest(tmp_path):
    rows = run(tmp_path, [r("R1", "S1", "P1", "SINGLE", bases=9_000_000_000), r("R2", "S1", "P1", bases=1_000_000_000),
                          r("R3", "S1", "P1", bases=2_000_000_000, tax="330879"), r("R4", "S2", "P1"), r("R5", "S9", "P1", tax="5062")])
    by = {x["sample_id"]: x for x in rows}
    assert set(by) == {"S1", "S2"}                       # other species left out; strain taxon rolled up
    assert by["S1"]["run_accession"] == "R3"             # paired beats single; most bases among paired
    # 30x of 30 Mb = 900 Mb at 300 bases per spot (S2: 3 Gb / 10 M spots) = 3,000,000 spots
    assert by["S2"]["target_spots"] == "3000000" and by["S2"]["full_coverage"] == "100.0"


def test_low_coverage_run_is_read_whole(tmp_path):
    rows = run(tmp_path, [r("R1", "S1", "P1", reads=1_000_000, bases=300_000_000)])
    assert rows[0]["target_spots"] == "1000000"


def test_flagged_study_dropped_and_local_rows(tmp_path):
    runs = [r(f"C{i}", f"SC{i}", "PX", title="Meiosis: progeny of a sexual cross", strain=f"progeny {i}") for i in range(4)]
    runs.append(r("R1", "S1", "P1", strain="AF293"))
    local = "AF293,SRA/x_[12].fastq.gz,Illumina\nNEW1,AF100/y_R?_001.fastq.gz,Illumina\nASM1,Pre_assembled/a.fasta,assembly\n"
    rows = run(tmp_path, runs, local)
    ids = [x["sample_id"] for x in rows]
    assert "SC0" not in ids and "S1" in ids
    loc = {x["sample_id"]: x for x in rows if x["source"] == "local"}
    assert set(loc) == {"local__AF293", "local__NEW1"}     # pre-assembled sample skipped
    assert loc["local__AF293"]["also_public"] == "S1" and loc["local__NEW1"]["layout"] == "PAIRED"
