"""datasets/reads/sra_survey/pilot/study_titles.py: titles for the studies of the wave folders and of the sketch-store reports."""
import csv
import gzip
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "datasets" / "reads" / "sra_survey" / "pilot" / "study_titles.py"


def tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(rows)


def test_wave_and_sketch_store_studies_get_titles(tmp_path):
    reads = tmp_path / "datasets" / "reads"
    pilot = reads / "sra_survey" / "pilot"
    tsv(pilot / "wave1" / "results_all_final.tsv", ["run", "study"], [["SRR1", "PRJNA1"], ["SRR2", ""]])
    tsv(reads / "sketch_store" / "reports" / "sp" / "runs_public.tsv", ["biosample", "run", "study"], [["SAMN9", "SRR9", "PRJNA9"], ["SAMN8", "SRR8", "PRJNA8"]])
    runs = tmp_path / "runs.tsv.gz"
    with gzip.open(runs, "wt") as fh:
        fh.write("run_accession\tstudy_accession\tstudy_title\nSRR1\tPRJNA1\tWave  study\nSRR9\tPRJNA9\tStore study\nSRR7\tPRJNA7\tNot used\n")
    out = tmp_path / "titles.tsv"
    subprocess.run([sys.executable, str(SCRIPT), str(pilot), str(runs), str(out)], check=True, capture_output=True)
    got = {r["study"]: r["title"] for r in csv.DictReader(open(out), delimiter="\t")}
    assert got == {"PRJNA1": "Wave study", "PRJNA8": "", "PRJNA9": "Store study"}
