#!/usr/bin/env python3
"""Run lists that top up typed species from a reads-availability snapshot: up to --cap runs per species in all (typed + new), at most --per-study new runs per study.

    expand_runs.py --snapshot snapshots/2026-11-06 --waves-root RESULTS/sra_reads_typing/waves --out OUT [--cap 300] [--per-study 25] [--actions expand unvalidated]

For each species whose action is in --actions, the untyped eligible runs of `untyped_runs.tsv.gz` (one run per isolate, already the run nearest 3 Gb) inside the 0.6-12 Gb
window are taken round-robin across studies (each study nearest 3 Gb first, then by accession), so no study supplies more than --per-study of them; the species stops at
--cap minus the isolates already typed. The panel columns (panel_dir, panel_prefix, idiomorphs) are copied from the species' own full run list in its wave
(WAVES_ROOT/<wave>/plan/<group>/runs_full.tsv), so the new runs are typed with exactly the panel that passed the gate.

Writes OUT/runs_<wave>.tsv (the run-list columns of ../sra_survey/pilot/CONTRACT.md; submit with `wave_driver.py submit-list --wave <wave> --runs ... --stage-name expand300`)
and OUT/expand_summary.tsv (per species: typed, available, selected, and what the window and the study cap left out). Standard library only.
"""
import argparse, collections, csv, gzip, os, sys

COLS = ["run_accession", "read_count", "group", "library_layout", "strain", "study_accession", "panel_dir", "panel_prefix", "idiomorphs"]


def read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pick(rows, need, per_study):
    """Round-robin across studies, nearest 3 Gb first inside a study; at most per_study from one study; at most need in all."""
    by = collections.defaultdict(list)
    for r in rows:
        by[r["study_accession"]].append(r)
    for q in by.values():
        q.sort(key=lambda r: (abs(int(r["base_count"] or 0) - 3e9), r["run_accession"]))
    order = sorted(by, key=lambda s: (-len(by[s]), s))       # the largest studies first in each round, so a cut at `need` falls on them
    out, rnd = [], 0
    while len(out) < need and rnd < per_study:
        took = False
        for s in order:
            if rnd < len(by[s]) and len(out) < need:
                out.append(by[s][rnd]); took = True
        if not took:
            break
        rnd += 1
    return out


def read_species_file(path):
    """Species names, one per line (the first tab-separated field); blank lines and lines starting with # are skipped. None when no file is given."""
    if not path:
        return None
    with open(path) as fh:
        return {ln.split("\t")[0].strip() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")}


def select_untyped(rows, actions, species=None):
    """{species: untyped rows} for the rows whose action is in `actions`, limited to `species` (a set) when given."""
    out = collections.defaultdict(list)
    for r in rows:
        if r["action"] in actions and (species is None or r["species"] in species):
            out[r["species"]].append(r)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--snapshot", required=True, help="a snapshot folder with species_availability.tsv and untyped_runs.tsv.gz")
    ap.add_argument("--waves-root", required=True, help="folder holding the wave directories (wave1, wave2, ...)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cap", type=int, default=300, help="most runs per species in all, typed ones included")
    ap.add_argument("--per-study", type=int, default=25, help="most NEW runs from one study (a guard against clonal outbreak sets)")
    ap.add_argument("--actions", nargs="+", default=["expand", "unvalidated"])
    ap.add_argument("--species-file", help="only these species (one name per line, # comments); default: every species with one of --actions")
    a = ap.parse_args(argv)

    sp = {r["species"]: r for r in read_tsv(os.path.join(a.snapshot, "species_availability.tsv"))}
    only = read_species_file(a.species_file)
    untyped = select_untyped(read_tsv(os.path.join(a.snapshot, "untyped_runs.tsv.gz")), a.actions, only)
    if only:
        for name in sorted(only - set(untyped)):
            print(f"{name}: no untyped run with an action in {a.actions} (gate {sp.get(name, {}).get('gate', 'not in the snapshot')})", file=sys.stderr)
    os.makedirs(a.out, exist_ok=True)
    lists, summary, missing = collections.defaultdict(list), [], []
    for name in sorted(untyped):
        rows = untyped[name]; wave = rows[0].get("wave", ""); group = rows[0]["group"]
        typed = int(sp.get(name, {}).get("typed_public") or 0)
        full = os.path.join(a.waves_root, wave, "plan", group, "runs_full.tsv")
        if not wave or not os.path.exists(full):
            missing.append((name, wave, group)); continue
        ref = read_tsv(full)
        panel = collections.Counter((r.get("panel_dir", ""), r.get("panel_prefix", ""), r.get("idiomorphs", "")) for r in ref).most_common(1)[0][0]
        inwin = [r for r in rows if r["in_select_window"] == "yes"]
        need = max(0, a.cap - typed)
        chosen = pick(inwin, need, a.per_study)
        uncapped = pick(inwin, need, 10 ** 9)
        for r in chosen:
            lists[wave].append(dict(run_accession=r["run_accession"], read_count=r["read_count"], group=group, library_layout=r["library_layout"], strain=r["strain"],
                                    study_accession=r["study_accession"], panel_dir=panel[0], panel_prefix=panel[1], idiomorphs=panel[2]))
        summary.append(dict(species=name, wave=wave, group=group, action=rows[0]["action"], gate=sp.get(name, {}).get("gate", ""), typed=typed, untyped=len(rows),
                            in_window=len(inwin), selected=len(chosen), left_by_study_cap=len(uncapped) - len(chosen), total_after=typed + len(chosen)))
    for wave, rows in lists.items():
        with open(os.path.join(a.out, f"runs_{wave}.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    with open(os.path.join(a.out, "expand_summary.tsv"), "w", newline="") as fh:
        cols = ["species", "wave", "group", "action", "gate", "typed", "untyped", "in_window", "selected", "left_by_study_cap", "total_after"]
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(summary)
    for name, wave, group in missing:
        print(f"no run list in a wave for {name} (wave {wave or '-'}, group {group}): skipped", file=sys.stderr)
    print(f"{sum(len(v) for v in lists.values())} runs for {len(summary)} species: " + ", ".join(f"{w} {len(v)}" for w, v in sorted(lists.items())))


if __name__ == "__main__":
    main()
