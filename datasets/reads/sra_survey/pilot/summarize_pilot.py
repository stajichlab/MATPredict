#!/usr/bin/env python3
"""Summarize one pilot or production array of stream_type_sra.slurm, per group (species): calls, agreement with an independent call,
failures, flags and timings. Reads only the small per-run files, so it can run any time and again (idempotent), and as a Slurm job after the array.

    summarize_pilot.py --runs RUNS.tsv --results RESULTS_DIR --out OUT_DIR

RUNS.tsv: the run list given to the array (columns run_accession, group; optional strain, truth_detect, truth_assembly, panel_prefix).
RESULTS_DIR: where the array wrote <run>.<panel_prefix>.tsv and <run>.run.json. Writes in OUT_DIR:
  results_all.tsv      one row per run: the run list, the call and its numbers, the timings, the state, agreement with the truth
  summary_by_group.tsv one row per group (species)
  discordant.tsv       typed runs whose call differs from the truth; failures.tsv: runs with no result and the state in their run.json
  SUMMARY.md           the same in words, with the numbers a reader needs to judge the run
"""
import argparse, collections, csv, glob, json, os, statistics


def read_tsv(path):
    with open(path, newline="") as fh:          # reads-type writes CRLF; a stray CR inside a line (a column appended after it) is dropped
        return list(csv.DictReader(fh.read().replace("\r", "").splitlines(), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", required=True); ap.add_argument("--results", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for r in read_tsv(a.runs):
        acc = r["run_accession"]; pfx = r.get("panel_prefix") or ""
        js = os.path.join(a.results, f"{acc}.run.json")
        j = json.load(open(js)) if os.path.exists(js) else {}
        tsvs = glob.glob(os.path.join(a.results, f"{acc}.{pfx or '*'}.tsv"))
        res = read_tsv(tsvs[0])[0] if tsvs else {}
        truth = (r.get("truth_detect") or "").strip()
        call = res.get("call", "")
        if not j:
            state = "not_run"
        else:
            state = j.get("state", "")
        typed = call in ("MAT1-1", "MAT1-2", "MATa", "MATalpha") or (call and call not in ("none", "low_depth", "both"))
        if truth and "," not in truth and call:
            verdict = "agree" if call == truth else ("both_vs_one" if call == "both" else ("untyped" if call in ("none", "low_depth") else "DISAGREE"))
        else:
            verdict = ""
        flags = res.get("flags", "")
        rows.append(dict(
            run=acc, group=r.get("group", ""), strain=r.get("strain", ""), study=r.get("study_accession", ""), layout=r.get("library_layout", j.get("layout", "")),
            spots=j.get("spots_total", ""), state=state, call=call, truth=truth, verdict=verdict,
            depth_1=res.get("MAT1-1_depth", ""), depth_2=res.get("MAT1-2_depth", ""), shared_depth=res.get("shared_depth", ""),
            breadth_1=res.get("MAT1-1_breadth", ""), breadth_2=res.get("MAT1-2_breadth", ""), flags=flags,
            sra_mb=round(j.get("sra_bytes", 0) / 1e6) if j else "", secs_download=j.get("secs_download", ""), secs_convert=j.get("secs_convert", ""),
            secs_type=j.get("secs_type", ""), truth_assembly=r.get("truth_assembly", "")))
    cols = list(rows[0].keys())
    with open(os.path.join(a.out, "results_all.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r["group"]].append(r)
    summ = []
    for g, rs in sorted(groups.items()):
        done = [r for r in rs if r["state"] == "done"]
        calls = collections.Counter(r["call"] for r in done)
        withT = [r for r in done if r["truth"] and "," not in r["truth"]]
        typedT = [r for r in withT if r["call"] in ("MAT1-1", "MAT1-2", "MATa", "MATalpha")]
        agree = [r for r in typedT if r["call"] == r["truth"]]
        mb = sum(r["sra_mb"] for r in done if r["sra_mb"] != ""); dl = sum(r["secs_download"] for r in done if r["secs_download"] != "")
        tasks = [r["secs_download"] + r["secs_convert"] + r["secs_type"] for r in done if r["secs_download"] != ""]
        summ.append(dict(group=g, runs=len(rs), done=len(done), failed=sum(1 for r in rs if r["state"].startswith("failed")), not_run=sum(1 for r in rs if r["state"] == "not_run"),
                         MAT1_1=calls.get("MAT1-1", 0), MAT1_2=calls.get("MAT1-2", 0), both=calls.get("both", 0), none=calls.get("none", 0), low_depth=calls.get("low_depth", 0),
                         typed_rate=round((calls.get("MAT1-1", 0) + calls.get("MAT1-2", 0)) / len(done), 3) if done else "",
                         with_truth=len(withT), typed_with_truth=len(typedT), agree=len(agree), disagree=len(typedT) - len(agree),
                         agreement=round(len(agree) / len(typedT), 3) if typedT else "", untyped_with_truth=len(withT) - len(typedT),
                         with_trace_flag=sum(1 for r in done if "trace" in r["flags"]), sra_gb=round(mb / 1000, 1),
                         mean_mb_per_s=round(mb / dl, 1) if dl else "", median_task_secs=statistics.median(tasks) if tasks else ""))
    with open(os.path.join(a.out, "summary_by_group.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(summ[0].keys()), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(summ)
    bad = [r for r in rows if r["verdict"] == "DISAGREE"]
    with open(os.path.join(a.out, "discordant.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(bad)
    fail = [r for r in rows if r["state"] != "done"]
    with open(os.path.join(a.out, "failures.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(fail)
    md = ["# Summary", "", f"Runs listed {len(rows)}; done {sum(1 for r in rows if r['state'] == 'done')}; failed {sum(1 for r in rows if r['state'].startswith('failed'))}; not run {sum(1 for r in rows if r['state'] == 'not_run')}.", "",
          "| group | runs | done | MAT1-1 | MAT1-2 | both | none | low_depth | with truth | typed | agree | disagree | untyped | GB | MB/s |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for s in summ:
        md.append(f"| {s['group']} | {s['runs']} | {s['done']} | {s['MAT1_1']} | {s['MAT1_2']} | {s['both']} | {s['none']} | {s['low_depth']} | {s['with_truth']} | {s['typed_with_truth']} | {s['agree']} | {s['disagree']} | {s['untyped_with_truth']} | {s['sra_gb']} | {s['mean_mb_per_s']} |")
    md += ["", "Truth = the `detect` idiomorph of a BFD assembly with the same BioSample (or a mating type in the strain name). `typed` = a one-idiomorph call; `untyped` = none or low_depth on a run that has a truth.",
           "Discordant runs: `discordant.tsv`; runs without a result: `failures.tsv`."]
    open(os.path.join(a.out, "SUMMARY.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
