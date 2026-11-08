#!/usr/bin/env python3
"""Round 3 merge (2026-11-08): fold the cap-700 stages, the Ophidiomyces panel fix and wave 5 into the committed wave tables; regate with the curator rulings.

Runs on the HPCC. Reads the committed tables from a checkout (PILOT) and the stage summaries from the waves folder (WAVES); writes new
results_all_final.tsv / gates_final.tsv (+ screen_gate_changes.tsv) per wave under OUT/<wave>/ for review before they are committed.

- expand700 + expand700_rerun: one row per run (done beats other states; the rerun wins among equals), merged with merge_screen.merge (new runs appended).
- Ophidiomyces (wave1): the wave6_ophidiomyces retype (panel with the 20 kb span) replaces the call where it typed; where it did not and the old
  panel did, the old call stays with rescue=old_panel.
- wave5_asco / wave5_mucoro: stages test + full combined; R. delemar untyped runs rescued by the Rhizopus genus panel (rescue=genus_panel).
- Gates: merge_screen.regate_species (wave_driver.gate_row: ZYGOSITY, RELEASED).

    merge_round3.py PILOT WAVES OUT
"""
import csv, os, sys
PILOT, WAVES, OUT = sys.argv[1:4]
sys.path.insert(0, PILOT)
import merge_screen as ms
import wave_driver as wd

COLS = ("call", "depth_1", "depth_2", "shared_depth", "breadth_1", "breadth_2", "flags")
rd = lambda p: list(csv.DictReader(open(p, newline=""), delimiter="\t")) if os.path.exists(p) else []
rank = lambda st: 2 if st == "done" else 1 if st in ("unavailable_at_ncbi", "failed_type", "done_untyped") else 0
typed = lambda r: r.get("state") == "done" and r.get("call") not in ("none", "low_depth", "", None)


def combine(paths):
    best = {}
    for p in paths:
        for r in rd(p):
            o = best.get(r["run"])
            if o is None or rank(r["state"]) >= rank(o["state"]): best[r["run"]] = r
    return list(best.values())


def write(wave, rows, gates, changes, species=None):
    d = os.path.join(OUT, wave); os.makedirs(d, exist_ok=True)
    cols = []
    for r in rows:
        for c in r:
            if c not in cols: cols.append(c)
    for name, data, cc in (("results_all_final.tsv", rows, cols), ("gates_final.tsv", gates, ["slug", "species", "decision", "reasons"]),
                           ("screen_gate_changes.tsv", changes, ["slug", "species", "previous", "decision", "reasons"])):
        with open(os.path.join(d, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, cc, delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(data)
    if species:
        with open(os.path.join(d, "species.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, ["species", "slug"], delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(species)
    from collections import Counter
    print(f"{wave}: {len(rows)} rows; gates {dict(Counter(g['decision'].split(':')[0] for g in gates))}; changed: " +
          ("; ".join(f"{c['species']} {c['previous'].split(':')[0]}->{c['decision'].split(':')[0]}" for c in changes) or "none"))


for wave in ("wave1", "wave2", "wave3_asco", "wave3_basidio", "wave3_gattii", "wave3_mucoro", "wave3_poryzae", "wave4", "wave4_groups"):
    final = rd(os.path.join(PILOT, wave, "results_all_final.tsv"))
    stage = combine([os.path.join(WAVES, wave, s, "summary", "results_all.tsv") for s in ("expand700", "expand700_rerun")])
    rows, added = ms.merge(final, stage) if stage else (final, 0)
    if wave == "wave1":
        new = {r["run"]: r for r in combine([os.path.join(WAVES, "wave6_ophidiomyces", s, "summary", "results_all.tsv") for s in ("retype", "retype700")])}
        n_new = n_old = 0
        for r in rows:
            x = new.get(r["run"])
            if r["group"] != "ophidiomyces_ophidiicola" or not x or x["state"] != "done": continue
            if typed(x):
                for k in COLS: r[k] = x.get(k, r.get(k, ""))
                r["rescue"] = ""; n_new += 1
            elif typed(r):
                r["rescue"] = "old_panel"; n_old += 1
            else:
                for k in COLS: r[k] = x.get(k, r.get(k, ""))
        print(f"  Ophidiomyces: {n_new} calls from the new panel, {n_old} kept from the old panel")
    gates, changes = ms.regate_species(rows, rd(os.path.join(PILOT, wave, "gates_final.tsv")))
    write(wave, rows, gates, changes)

# wave 5: test + full; R. delemar genus rescue
resc = {r["run"]: r for r in rd(os.path.join(WAVES, "genus_rescue2", "rdelemar", "summary", "results_all.tsv"))}
for wave in ("wave5_asco", "wave5_basidio", "wave5_mucoro"):
    rows = combine([os.path.join(WAVES, wave, s, "summary", "results_all.tsv") for s in ("test", "full")])
    for r in rows:
        x = resc.get(r["run"])
        if x and r["state"] == "done" and r["call"] in ("none", "low_depth") and typed(x):
            for k in COLS: r[k] = x.get(k, r.get(k, ""))
            r["rescue"] = "genus_panel"
    gates = rd(os.path.join(WAVES, wave, "gates.tsv"))
    by = {}
    for r in rows: by.setdefault(r["group"], []).append(r)
    out_g = [wd.gate_row(g["slug"], g["species"], by[g["slug"]]) if by.get(g["slug"]) else g for g in gates]
    write(wave, rows, out_g, [], species=rd(os.path.join(WAVES, wave, "species.tsv")))
