"""datasets/reads/sra_survey/pilot/cause_evidence.py: assembly evidence per species (one idiomorph, both in one locus, both at separate loci) beside the gate decision."""
import csv, subprocess, sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent.parent / "datasets" / "reads" / "sra_survey" / "pilot" / "cause_evidence.py"


def tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(rows)


def test_genomes_are_counted_by_how_their_mat_loci_look(tmp_path):
    pilot = tmp_path / "pilot"; rel = tmp_path / "rel"
    tsv(pilot / "wave1" / "gates_final.tsv", ["slug", "species", "decision", "reasons"],
        [["sp_h", "Sp homo", "PROVISIONAL_HOMOTHALLIC", "runs=4 done=4 typed=0 both=4 none/low=0"], ["sp_ok", "Sp ok", "PASS", "runs=4"], ["sp_n", "Sp named", "UNVALIDATED:x", "runs=1"]])
    tsv(pilot / "causes.tsv", ["species", "cause", "build_cause", "basis"], [["Sp named", "silent_cassettes", "", "b"]])
    tsv(rel / "genomes.tsv", ["genome", "species"], [["g1", "Sp homo"], ["g2", "Sp homo"], ["g3", "Sp homo"], ["g4", "Sp homo"], ["g5", "Sp ok"], ["g6", "Sp named"]])
    tsv(rel / "loci.tsv", ["genome", "idiomorph", "locus_class"],
        [["g1", "MAT1-1", "mat_locus"], ["g1", "MAT1-2", "mat_locus"],      # both, separate loci
         ["g2", "undetermined", "homothallic_candidate"],                     # both in one locus
         ["g3", "MAT1-2", "mat_locus"],                                       # one idiomorph; g4 has no locus
         ["g5", "MAT1-1", "mat_locus"], ["g6", "MATa", "mat_locus"]])
    out = tmp_path / "ev.tsv"
    subprocess.run([sys.executable, str(SCRIPT), "--pilot", str(pilot), "--release", str(rel), "--out", str(out)], check=True, capture_output=True)
    rows = {r["species"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert set(rows) == {"Sp homo", "Sp named"}                               # the PASS species is not a target; a species in causes.tsv is
    h = rows["Sp homo"]
    assert (h["genomes"], h["genomes_with_mat"], h["single"], h["both_one_locus"], h["both_separate"]) == ("4", "3", "1", "1", "1")
    assert h["reads"] == "done 4, single 0, both 4" and h["labels"].startswith("MAT1-2: 2")
