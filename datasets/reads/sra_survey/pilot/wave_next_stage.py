#!/usr/bin/env python3
"""After the rescue passes: which species pass the final gate but did not pass the test-stage gate, and which of their runs (up to --cap per species, the test runs first) still have no result? Writes the run list to --out.

    wave_next_stage.py --wave WAVE_DIR [--cap 50] [--out full3_runs.tsv]      then:  wave_driver.py submit-list --wave WAVE_DIR --runs full3_runs.tsv --stage-name full3 --env FALLBACK_K=21 --env SKETCH=1
"""
import argparse, csv, os

OK = ("PASS", "PROVISIONAL", "PROVISIONAL_HOMOTHALLIC")
rd = lambda p: list(csv.DictReader(open(p, newline="").read().replace("\r", "").splitlines(), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--wave", required=True); ap.add_argument("--cap", type=int, default=50); ap.add_argument("--out")
    ap.add_argument("--final-gates", default="gates_final.tsv"); a = ap.parse_args(); W = a.wave
    old = {r["slug"]: r["decision"] for r in rd(W + "/gates.tsv")}; fin = {r["slug"]: r["decision"] for r in rd(os.path.join(W, a.final_gates))}
    new = [s for s, d in fin.items() if d in OK and old.get(s) not in OK]
    done = {f.split(".")[0] for f in os.listdir(W + "/out") if f.endswith(".run.json")}
    todo = []
    for s in new:
        rows = rd(f"{W}/plan/{s}/runs_full.tsv")[:a.cap]; left = [r for r in rows if r["run_accession"] not in done]; todo += left
        print(s, fin[s], "runs", len(rows), "to do", len(left))
    print("passing at the final gate:", sum(1 for d in fin.values() if d in OK), "| at the test gate:", sum(1 for d in old.values() if d in OK), "| newly passing:", len(new), "| runs to do:", len(todo))
    if a.out and todo:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, list(todo[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(todo)


if __name__ == "__main__":
    main()
