#!/usr/bin/env python3
"""Add the contamination-screen columns of a screened full stage to a wave's final results table (2026-11-04).

A run of the stage keeps its existing row, including a rescued call (`rescue`), and gains the screen columns; a run that the final table does not have is appended as it came
out of the stage. Rows of other runs get blank screen columns. Standard library only.

    merge_screen.py --final results_all_final.tsv --stage full_unvalidated/summary/results_all.tsv --out results_all_final.tsv
"""
import argparse, collections, csv, json, os, statistics, sys

SCREEN = ("screen", "screen_genus", "screen_second", "screen_hashes", "screen_explained", "screen_deeper")


def rd(p):
    return list(csv.DictReader(open(p, newline="").read().replace("\r", "").splitlines(), delimiter="\t"))


def merge(final, stage, species_of=lambda run: ""):
    """(rows, added): `final` rows with screen columns, then the stage runs the final table lacked. A stage run screened `no_call` under the first rule
    (no genus call) is `ok` when GmG gave a species call, as screen_reads.py now has it; `species_of(run)` gives that call."""
    stage = [dict(s, screen="ok") if s.get("screen") == "no_call" and species_of(s["run"]) else s for s in stage]
    by = {r["run"]: r for r in stage}
    rows, seen = [], set()
    for r in final:
        r = dict(r)
        s = by.get(r["run"], {})
        for k in SCREEN:
            r[k] = s.get(k, "")
        seen.add(r["run"])
        rows.append(r)
    added = [dict(s, rescue="") for s in stage if s["run"] not in seen]
    return rows + added, len(added)


def regate_species(rows, gates, decide=None):
    """Decide again, from the merged table, the species that gained runs. Returns (gates, changes): the gates with the new numbers (`reasons`) and decision of those
    species, and one row per species whose decision changed. `decide` is wave_driver.decide."""
    by = collections.defaultdict(list)
    for r in rows:
        by[r["group"]].append(r)
    changes = []
    for g in gates:
        rs = by.get(g["slug"], [])
        if not rs or not g["reasons"].startswith("runs="):
            continue
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from wave_driver import gate_row          # applies the curator rulings (ZYGOSITY counting, RELEASED species) as the gate does
        row = gate_row(g["slug"], g["species"], rs)
        new, reasons = row["decision"], row["reasons"]
        if reasons == g["reasons"] and new == g["decision"]:
            continue
        if new.split(":")[0] != g["decision"].split(":")[0]:
            changes.append({"slug": g["slug"], "species": g["species"], "previous": g["decision"], "decision": new, "reasons": reasons})
        g["decision"], g["reasons"] = new, reasons
    return gates, changes


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--final", required=True); ap.add_argument("--stage", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--gates", help="the wave's gates_final.tsv: decide again the species that gained runs (rewritten in place; changes go to screen_gate_changes.tsv beside it)")
    ap.add_argument("--screen-dir", help="folder of <run>.screen.json: re-reads the species call of a no_call run")
    a = ap.parse_args(argv)
    final = rd(a.final)
    cols = list(final[0].keys())
    cols += [k for k in SCREEN if k not in cols]
    def species_of(run):
        try:
            return json.load(open(os.path.join(a.screen_dir, f"{run}.screen.json"))).get("species") or ""
        except (OSError, TypeError, ValueError):
            return ""
    rows, added = merge(final, rd(a.stage), species_of)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    if a.gates:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        gates, changes = regate_species(rows, rd(a.gates))
        for path, data, cols in ((a.gates, gates, ["slug", "species", "decision", "reasons"]),
                                 (os.path.join(os.path.dirname(a.gates), "screen_gate_changes.tsv"), changes, ["slug", "species", "previous", "decision", "reasons"])):
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(data)
        print(f"gates: {len(changes)} decisions changed: " + "; ".join(f"{c['species']} {c['previous'].split(':')[0]} -> {c['decision'].split(':')[0]}" for c in changes))
    print(f"{len(rows)} rows ({added} added), {sum(1 for r in rows if r['screen'])} screened")


if __name__ == "__main__":
    main()
