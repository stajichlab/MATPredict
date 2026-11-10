"""Fold the thin-evidence stage (`thin`) into the committed wave tables of a repo checkout (curator 2026-11-09).

For each wave: the stage's runs are appended to results_all_final.tsv (runs already there keep their row unchanged), and the gates_final.tsv
rows of the re-planned species are replaced by their rows of the wave's gates_thin.tsv; changes go to thin_changes.tsv. A wave with no committed folder
(wave5_basidio) gets one: species.tsv, gates_final.tsv (all species of gates_thin.tsv), results_all_final.tsv (the stage's runs).
    merge_thin.py REPO
"""
import csv, os, sys
REPO = sys.argv[1]
PILOT = f"{REPO}/datasets/reads/sra_survey/pilot"
W = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves"
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-09_thin_evidence_replan"
sys.path.insert(0, PILOT)
import merge_screen as ms


def rd(p):
    return ms.rd(p) if os.path.exists(p) else []


def wr(p, rows, cols):
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)


GCOLS = ["slug", "species", "decision", "reasons"]
for w in ("wave1", "wave2", "wave4", "wave5_asco", "wave5_basidio", "wave5_mucoro"):
    thin_slugs = {r["group"] for r in rd(f"{R}/{w}_full.tsv")}
    stage = rd(f"{W}/{w}/thin/summary/results_all.tsv")
    gthin = {r["slug"]: r for r in rd(f"{W}/{w}/gates_thin.tsv")}
    d = f"{PILOT}/{w}"; new_folder = not os.path.isdir(d); os.makedirs(d, exist_ok=True)
    final = rd(f"{d}/results_all_final.tsv")
    if final:      # append only: existing rows keep every column (merge_screen.merge would blank the screen columns of runs outside this stage)
        cols = list(final[0].keys()); seen = {r["run"] for r in final}
        new = [dict(s, rescue="") for s in stage if s["run"] not in seen]
        rows, added = final + new, len(new)
    else:
        cols = list(stage[0].keys()) + (["rescue"] if "rescue" not in stage[0] else []); rows, added = [dict(s, rescue="") for s in stage], len(stage)
    wr(f"{d}/results_all_final.tsv", rows, cols)
    gates = rd(f"{d}/gates_final.tsv"); changes = []
    if new_folder or not gates:
        gates = list(gthin.values())
        wr(f"{d}/species.tsv", [{"species": g["species"], "slug": g["slug"]} for g in gates], ["species", "slug"])
        changes = [dict(slug=g["slug"], species=g["species"], previous="(new folder)", decision=g["decision"], reasons=g["reasons"]) for g in gates if g["slug"] in thin_slugs]
    else:
        for i, g in enumerate(gates):
            if g["slug"] in thin_slugs and g["slug"] in gthin:
                n = gthin[g["slug"]]; changes.append(dict(slug=g["slug"], species=g["species"], previous=g["decision"], decision=n["decision"], reasons=n["reasons"])); gates[i] = {k: n[k] for k in GCOLS}
    wr(f"{d}/gates_final.tsv", gates, GCOLS)
    wr(f"{d}/thin_changes.tsv", changes, ["slug", "species", "previous", "decision", "reasons"])
    print(w, f"runs added {added}", f"gates changed {len(changes)}", "(new folder)" if new_folder else "")
    for c in changes: print("  ", c["species"], c["previous"], "->", c["decision"])
