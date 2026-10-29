#!/usr/bin/env python3
"""Run raw-read MAT typing for a WAVE of species: plan, test, gate, full run, report. Idempotent; every step can be repeated; state is files.

A wave directory (WAVE) holds:
  species.tsv                       input: species<TAB>slug (one row per species)
  plan/<slug>/                      panel/ (lineage panels L*/), runs_test.tsv, runs_full.tsv, qc_panel.json, status (planned | held:<reason>)
  out/                              ONE folder of per-run results for all species (<run>.<prefix>.tsv, <run>.run.json); a finished run is never redone
  test/, full/                      run lists, job ids and summaries of the two stages
  gates.tsv                         the decision per species and the numbers behind it
  WAVE_REPORT.md                    across species

Commands (all take --wave DIR and the config options below):
  plan-one SPECIES     build the lineage panels, check them (panel gate), choose test and full run lists for ONE species
  submit-plan          queue `plan-one` for every species of species.tsv as an array (panels are built on compute nodes)
  submit-test          merge the test runs of the planned species and queue them (split over the exfab and stajichlab partitions)
  gate                 summarize the test results and decide per species: PASS, PROVISIONAL or HOLD:<reasons>   (gates below)
  submit-full          queue the remaining runs of the species that passed the gate
  report               write WAVE_REPORT.md from the summaries and the gates
  status               counts of species by state, and of runs done and failed

Gates. Panel (before any download): lineage panel L1 exists; at least 2 loci of each idiomorph went into the panel; core k-mers per idiomorph >= 300 and within a factor 10
of each other; shared flank k-mers >= 2000; at least 4 selectable runs (the survey's own minimum; a species with fewer than 8 runs is typed whole and can only be PROVISIONAL). Test (after 20 runs): at least 90% of runs finished (runs not at NCBI are not counted); typed
(one idiomorph) >= 75% of finished runs and `both` <= 20%; median shared-flank depth >= 3; with >= 5 typed runs that have a `detect` truth: agreement >= 95% and at most one
disagreement; with fewer than 5 such runs the species is PROVISIONAL (typed >= 85% required). A species with at least half of its runs `both` and typed+both >= 85% is PROVISIONAL_HOMOTHALLIC (`both` is the expected answer for a homothallic species). Thresholds are in GATE below.
"""
import argparse, collections, csv, json, os, re, statistics, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
P = "/bigdata/stajichlab/jstajich/projects/MATPredict"
PY = f"{P}/.pixi/envs/default/bin/python"
GATE = dict(min_loci=2, min_core_kmers=300, max_core_ratio=10, min_shared_kmers=2000, min_runs=4,
            min_done=0.90, min_typed=0.75, max_both=0.20, min_shared_depth=3.0, min_truth_runs=5, min_agree=0.95, max_disagree=1, provisional_typed=0.85, homothallic_both=0.5)
TSV = dict(delimiter="\t", lineterminator="\n")


def rd(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh.read().replace("\r", "").splitlines(), delimiter="\t"))


def wr(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, **TSV); w.writeheader(); w.writerows(rows)


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, **kw)


def slug(sp): return re.sub(r"[^a-z0-9]+", "_", sp.lower()).strip("_")


def cfg(a):
    return dict(wave=a.wave, inputs=a.inputs)


# ------------------------------------------------------------------------------------------------ plan
# Genera whose MAT loci are only ever reported as a gene (locus class idiomorph_gene_only) and are short: the Taphrinomycotina (fission yeast, Pneumocystis, Taphrina ...). They are planned with that class
# allowed and a lower minimum span (the loci of S. japonicus are 1 to 2.6 kb), and an idiomorph missing from the species is borrowed from the genus.
LENIENT_GENERA = ("Schizosaccharomyces", "Pneumocystis", "Taphrina", "Protomyces", "Saitoella", "Neolecta")


def plan_one(a):
    sp = a.species; sl = slug(sp); d = os.path.join(a.wave, "plan", sl); os.makedirs(d, exist_ok=True)
    I = a.inputs; prefix = sl.split("_")[0]; panel = os.path.join(d, "panel")
    log = open(os.path.join(d, "plan.log"), "w"); lenient = sp.split()[0] in LENIENT_GENERA
    def hold(reason):
        open(os.path.join(d, "status"), "w").write("held:" + reason + "\n"); print(sp, "HELD", reason); sys.exit(0)
    r = subprocess.run([PY, f"{HERE}/panel_from_assemblies.py", "--species", sp, "--genomes", f"{I}/asco_genomes.tsv", "--loci", f"{I}/asco_loci.tsv", "--reports", f"{I}/reports_all.tar.zst",
                        "--genome-dir", f"{I}/genomes", "--out", panel, "--prefix", prefix, "--min-cluster-size", "1", "--max-reps", "6", "--per-lineage", "--jobs", "4"] + (["--min-span", "1000", "--classes", "mat_locus", "partial_locus", "idiomorph_gene_only"] if lenient else ["--min-span", "4000"]) + (["--borrow", sp.split()[0] + " "] if getattr(a, "borrow", False) or lenient else []),
                       stdout=log, stderr=log)
    if r.returncode or not os.path.isdir(os.path.join(panel, "L1")): hold("panel_build_failed")
    sys.path.insert(0, f"{P}/.claude/worktrees/run-b317ece/src")
    from MATPredict.reads.panel import Panel
    nf = os.path.join(panel, f"{prefix}_idiomorphs.txt"); names = open(nf).read().split() if os.path.exists(nf) else ["MAT1-1", "MAT1-2"]
    qc = dict(species=sp, prefix=prefix, lineages={}, idiomorphs=names, borrowed=bool(getattr(a, "borrow", False)))
    lin = rd(os.path.join(panel, f"{prefix}_lineages.tsv")) if os.path.exists(os.path.join(panel, f"{prefix}_lineages.tsv")) else []
    qc["loci_per_idiomorph"] = {i: sum(int(x[i]) for x in lin) for i in names}
    for ld in sorted(os.listdir(panel)):
        if re.fullmatch(r"L\d+", ld):
            try:
                pn = Panel.from_fastas({i: f"{panel}/{ld}/{prefix}_{i}.fasta" for i in names}, k=31, min_run=100)
                qc["lineages"][ld] = dict(core={k: len(v) for k, v in pn.unique.items()}, shared=len(pn.shared))
            except Exception as ex:
                qc["lineages"][ld] = dict(error=str(ex)[:120])
    json.dump(qc, open(os.path.join(d, "qc_panel.json"), "w"), indent=1)
    l1 = qc["lineages"].get("L1", {})
    if "no unique k-mers" in l1.get("error", ""): hold("panel_no_unique_kmers")      # the lineage's loci carry both idiomorphs' sequence (homothallic or switching species): presence of an idiomorph cannot be typed from reads
    if "error" in l1 or not l1: hold("panel_L1_error")
    if min(qc["loci_per_idiomorph"].values()) < GATE["min_loci"]: hold("too_few_loci")
    c = list(l1["core"].values())
    if min(c) < GATE["min_core_kmers"]: hold("core_kmers_low")
    if max(c) / max(1, min(c)) > GATE["max_core_ratio"]: hold("core_kmers_unbalanced")
    if l1["shared"] < GATE["min_shared_kmers"]: hold("shared_kmers_low")
    common = ["--species", sp, "--group", sl, "--panel-dir", panel, "--panel-prefix", prefix, "--runs", f"{I}/fungi_wgs_illumina.tsv", "--meta", f"{I}/fungi_meta.tsv",
              "--names", *names] + (["--classes", "mat_locus", "partial_locus", "idiomorph_gene_only"] if lenient else []) + ["--studies", f"{I}/studies.tsv", "--assembly-summary", f"{I}/asm_summary_fungi.txt", "--genomes", f"{I}/asco_genomes.tsv", "--loci", f"{I}/asco_loci.tsv"]
    r = subprocess.run([PY, f"{HERE}/select_runs.py", *common, "--n-truth", "12", "--n-free", "8", "--out", f"{d}/runs_test.tsv"], capture_output=True, text=True); log.write(r.stdout + r.stderr)
    if r.returncode or not os.path.exists(f"{d}/runs_test.tsv"): hold("run_selection_failed")
    test = rd(f"{d}/runs_test.tsv")
    if len(test) < GATE["min_runs"]: hold(f"too_few_runs_{len(test)}")
    r = subprocess.run([PY, f"{HERE}/select_runs.py", *common, "--n-truth", "25", "--n-free", "25", "--include", ",".join(x["run_accession"] for x in test), "--out", f"{d}/runs_full.tsv"], capture_output=True, text=True); log.write(r.stdout + r.stderr)
    open(os.path.join(d, "status"), "w").write("planned\n"); print(sp, "planned:", len(test), "test runs;", len(rd(f"{d}/runs_full.tsv")), "full runs")


# ------------------------------------------------------------------------------------------------ submit helpers
def sbatch(args, **kw):
    r = sh(["sbatch", "--parsable", *args]); 
    if r.returncode: sys.exit(f"sbatch failed: {r.stderr.strip()}")
    return r.stdout.strip().split(";")[0]


def submit_stage(a, stage, rows):
    """rows: run-list rows of ALL species for this stage; split over the two partitions; arrays of at most 2000 tasks; retry pass; summary job."""
    d = os.path.join(a.wave, stage); os.makedirs(d, exist_ok=True); os.makedirs(os.path.join(d, "logs"), exist_ok=True)
    out = a.out_dir or os.path.join(a.wave, "out"); os.makedirs(out, exist_ok=True)
    cols = list(rows[0].keys()); allp = os.path.join(d, "runs_all.tsv"); wr(allp, rows, cols)
    parts = [("exfab", "exfab"), ("stajichlab", "stajichlab")]
    ids = []; jobs = []
    for pi, (part, acct) in enumerate(parts):
        mine = rows[pi::2]
        for ci in range(0, len(mine), 2000):
            chunk = mine[ci:ci + 2000]; f = os.path.join(d, f"runs_{part}_{ci // 2000}.tsv"); wr(f, chunk, cols)
            for attempt in ("run", "retry"):
                dep = ["--dependency", f"afterany:{ids[-1]}"] if attempt == "retry" else []
                jid = sbatch(["-p", part, "-A", acct, "-c", str(a.cpus), "--mem", a.mem, "-J", f"wave-{stage}-{part}-{attempt}", "--array", f"1-{len(chunk)}%{a.conc}", *dep, "-o", f"{d}/logs/%x_%A_%a.out",
                              "--export", f"ALL,WORKER_DIR={HERE},TARGET_SPOTS=2000000,WINDOWS=8" + "".join("," + e for e in a.env), f"{HERE}/stream_type_sra.slurm", f, "none", "none", out]) if not a.dry_run else f"dry{len(ids)}"
                ids.append(jid); jobs.append((stage, part, attempt, jid, len(chunk)))
    last = ids[-1]
    sm = sbatch(["-p", "exfab", "-A", "exfab", "-c", "1", "--mem", "4gb", "-t", "0:30:00", "-J", f"wave-{stage}-summary", "--dependency", "afterany:" + ":".join(ids), "-o", f"{d}/logs/%x_%j.out",
                 "--wrap", f"{PY} {HERE}/summarize_pilot.py --runs {allp} --results {out} --out {d}/summary"]) if not a.dry_run else "drysummary"
    jobs.append((stage, "-", "summary", sm, 0))
    if not a.dry_run:
        with open(os.path.join(a.wave, "jobs.tsv"), "a") as fh:
            for j in jobs: fh.write("\t".join(map(str, j)) + "\n")
    print(f"{stage}: {len(rows)} runs, arrays {ids}, summary {sm}")


def submit_list(a):
    rows = rd(a.runs)
    if not rows: sys.exit("empty run list")
    submit_stage(a, a.stage_name, rows)


def planned(a):
    out = []
    for r in rd(os.path.join(a.wave, "species.tsv")):
        sf = os.path.join(a.wave, "plan", r["slug"], "status")
        if os.path.exists(sf) and open(sf).read().strip() == "planned": out.append(r["slug"])
    return out


def submit_plan(a):
    sp = rd(os.path.join(a.wave, "species.tsv")); os.makedirs(os.path.join(a.wave, "plan"), exist_ok=True)
    if a.retry_held or a.retry_status:      # plan again the species held for too few runs (after the gate changed), or whose status starts with --retry-status (for example held:panel_build_failed)
        pre = a.retry_status or "held:too_few_runs"
        sp = [r for r in sp if os.path.exists(os.path.join(a.wave, "plan", r["slug"], "status")) and open(os.path.join(a.wave, "plan", r["slug"], "status")).read().startswith(pre)]; os.makedirs(os.path.join(a.wave, "plan_logs"), exist_ok=True)
    lst = os.path.join(a.wave, f"plan_species_{int(__import__('time').time())}.txt"); open(lst, "w").write("".join(r["species"] + "\n" for r in sp))      # one list per submission: a pending array reads its own list, not a later one
    cmd = (f"export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/bin:$PATH; sp=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p {lst}); {PY} {HERE}/wave_driver.py plan-one --wave {a.wave} --inputs {a.inputs} --species \"$sp\"" + (" --borrow" if a.borrow else ""))
    jid = sbatch(["-p", "exfab", "-A", "exfab", "-c", "4", "--mem", "16gb", "-t", "2:00:00", "-J", "wave-plan", "--array", f"1-{len(sp)}%{a.conc}", "-o", f"{a.wave}/plan_logs/%x_%A_%a.out", "--wrap", cmd]) if not a.dry_run else "dryplan"
    print("plan array", jid, "for", len(sp), "species")


def submit_test(a):
    rows = []
    for sl in planned(a): rows += rd(os.path.join(a.wave, "plan", sl, "runs_test.tsv"))
    if not rows: sys.exit("no planned species")
    submit_stage(a, "test", rows)


# ------------------------------------------------------------------------------------------------ gate
def gate(a):
    res = rd(os.path.join(a.wave, a.stage, "summary", "results_all.tsv"))
    by = collections.defaultdict(list)
    for r in res: by[r["group"]].append(r)
    rows = []
    for r0 in rd(os.path.join(a.wave, "species.tsv")):
        sl = r0["slug"]; sf = os.path.join(a.wave, "plan", sl, "status"); st = open(sf).read().strip() if os.path.exists(sf) else "not_planned"
        if st != "planned":
            rows.append(dict(slug=sl, species=r0["species"], decision="HOLD:" + st.replace("held:", ""), reasons=st)); continue
        rs = by.get(sl, [])
        if not rs: rows.append(dict(slug=sl, species=r0["species"], decision="PENDING", reasons="no test results yet")); continue
        counted = [x for x in rs if x["state"] != "unavailable_at_ncbi"]
        done = [x for x in counted if x["state"] == "done"]; why = []
        typed = [x for x in done if x["call"] not in ("both", "none", "low_depth", "")]; both = [x for x in done if x["call"] == "both"]
        tr = [x for x in typed if x["truth"] and "," not in x["truth"]]; agree = sum(1 for x in tr if x["call"] == x["truth"])
        sd = [float(x["shared_depth"]) for x in done if x["shared_depth"] not in ("", None)]
        msd = statistics.median(sd) if sd else 0
        if not counted or len(done) / len(counted) < GATE["min_done"]: why.append("many_failed")
        tf = len(typed) / max(1, len(done)); bf = len(both) / max(1, len(done))
        homo = bf >= GATE["homothallic_both"] and (len(typed) + len(both)) / max(1, len(done)) >= GATE["provisional_typed"]    # a homothallic species: `both` is the expected answer
        if not homo:
            if tf < GATE["min_typed"]: why.append(f"typed_{tf:.2f}")
            if bf > GATE["max_both"]: why.append(f"both_{bf:.2f}")
        if msd < GATE["min_shared_depth"]: why.append(f"shared_depth_{msd:.1f}")
        prov = len(tr) < GATE["min_truth_runs"]
        if homo: prov = True
        if not prov:
            if agree / len(tr) < GATE["min_agree"] or len(tr) - agree > GATE["max_disagree"]: why.append(f"agree_{agree}of{len(tr)}")
        elif not homo and tf < GATE["provisional_typed"]: why.append(f"provisional_typed_{tf:.2f}")
        dec = "HOLD:" + "+".join(why) if why else ("PROVISIONAL_HOMOTHALLIC" if homo else "PROVISIONAL" if prov else "PASS")
        rows.append(dict(slug=sl, species=r0["species"], decision=dec, reasons=f"runs={len(rs)} done={len(done)} typed={len(typed)} both={len(both)} none/low={len(done)-len(typed)-len(both)} truth_typed={len(tr)} agree={agree} median_shared_depth={msd:.1f}"))
    wr(os.path.join(a.wave, "gates.tsv" if a.stage == "test" else f"gates_{a.stage}.tsv"), rows, ["slug", "species", "decision", "reasons"])
    print(collections.Counter(r["decision"].split(":")[0] for r in rows))


def submit_full(a):
    g = {r["slug"]: r["decision"] for r in rd(os.path.join(a.wave, "gates.tsv"))}
    rows = []
    for sl, dec in g.items():
        if dec in ("PASS", "PROVISIONAL", "PROVISIONAL_HOMOTHALLIC"): rows += rd(os.path.join(a.wave, "plan", sl, "runs_full.tsv"))[:a.cap]     # the test runs come first; at most --cap runs per species
    if not rows: sys.exit("no species passed the gate")
    submit_stage(a, "full", rows)


# ------------------------------------------------------------------------------------------------ report / status
def report(a):
    g = {r["slug"]: r for r in rd(os.path.join(a.wave, "gates_full.tsv" if os.path.exists(os.path.join(a.wave, "gates_full.tsv")) else "gates.tsv"))}
    stage = "full" if os.path.exists(os.path.join(a.wave, "full", "summary", "results_all.tsv")) else "test"
    res = rd(os.path.join(a.wave, stage, "summary", "results_all.tsv")); by = collections.defaultdict(list)
    for r in res: by[r["group"]].append(r)
    L = [f"# Wave report ({stage} stage)", "", f"Species planned {sum(1 for v in g.values() if v['decision'] not in ('HOLD:panel_build_failed',))}; decisions: {dict(collections.Counter(v['decision'].split(':')[0] for v in g.values()))}.", "",
         "| species | decision | runs | MAT1-1 | MAT1-2 | both | none | low_depth | unavailable | truth typed | agree | studies |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for sl, v in sorted(g.items()):
        rs = by.get(sl, []); c = collections.Counter(x["call"] for x in rs if x["state"] == "done"); n1, n2 = ((rs[0].get("idiomorphs") if rs else "") or "MAT1-1 MAT1-2").split(); c["MAT1-1"], c["MAT1-2"] = c[n1], c[n2]; un = sum(1 for x in rs if x["state"] == "unavailable_at_ncbi")
        tr = [x for x in rs if x["state"] == "done" and x["truth"] and "," not in x["truth"] and x["call"] not in ("both", "none", "low_depth", "")]
        L.append(f"| {v['species']} | {v['decision']} | {len(rs)} | {c['MAT1-1']} | {c['MAT1-2']} | {c['both']} | {c['none']} | {c['low_depth']} | {un} | {len(tr)} | {sum(1 for x in tr if x['call'] == x['truth'])} | {len({x['study'] for x in rs})} |")
    L += ["", "Calls are from a subsample of each run (2 million spots, 8 windows) typed against lineage panels cut from BFD assemblies; `truth` = the `detect` idiomorph of an assembly with the same BioSample.",
          "Counts of strains are NOT population frequencies: related or clonal strains are not removed yet (4 runs per study at most is the only guard)."]
    open(os.path.join(a.wave, "WAVE_REPORT.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L[:6]))


def status(a):
    c = collections.Counter()
    for r in rd(os.path.join(a.wave, "species.tsv")):
        sf = os.path.join(a.wave, "plan", r["slug"], "status")
        c[open(sf).read().strip() if os.path.exists(sf) else "not_planned"] += 1
    print("species:", dict(c))
    out = os.path.join(a.wave, "out")
    if os.path.isdir(out):
        def state(f):
            try: return json.load(open(os.path.join(out, f))).get("state")
            except ValueError: return "failed_unfinished"
        st = collections.Counter(state(f) for f in os.listdir(out) if f.endswith(".run.json"))
        print("runs:", dict(st))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cmd", choices=["plan-one", "submit-plan", "submit-test", "gate", "submit-full", "report", "status", "submit-list"])
    ap.add_argument("--wave", required=True); ap.add_argument("--inputs", default=f"{P}/results/sra_reads_typing/wave_inputs")
    ap.add_argument("--species"); ap.add_argument("--conc", type=int, default=60, help="tasks running at once per array (the curator lifted the CPU limit for this work, 2026-10-27)")
    ap.add_argument("--cpus", type=int, default=1, help="CPUs per run task: every step is single-threaded (measured 2026-10-27), and the exfab account is limited to 32 CPUs, so 1 CPU gives 32 tasks at once instead of 8")
    ap.add_argument("--mem", default="6gb", help="memory per run task (peak measured 4.9 GB)")
    ap.add_argument("--cap", type=int, default=50, help="most runs per species in the full stage")
    ap.add_argument("--stage", default="test", help="gate: which stage's summary to judge (a folder under the wave directory with summary/results_all.tsv) (the full stage has more runs per species)")
    ap.add_argument("--env", action="append", default=[], help="extra KEY=VALUE for the run tasks (for example READS_K=21); repeatable")
    ap.add_argument("--runs", help="submit-list: a run-list TSV with panel_dir and panel_prefix columns")
    ap.add_argument("--stage-name", default="list", help="submit-list: stage folder name under the wave directory")
    ap.add_argument("--out-dir", help="submit-list: results folder (default <wave>/out; the worker skips a run whose result is already there, so a different analysis of the same runs, for example sketches, needs its own)")
    ap.add_argument("--retry-held", action="store_true", help="submit-plan: only the species held for too few runs")
    ap.add_argument("--retry-status", help="submit-plan: only the species whose plan status starts with this (for example held:panel_build_failed)")
    ap.add_argument("--borrow", action="store_true", help="plan-one, submit-plan: an idiomorph with no locus in the species takes its loci from the other species of the genus")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    {"plan-one": plan_one, "submit-plan": submit_plan, "submit-test": submit_test, "gate": gate, "submit-full": submit_full, "report": report, "status": status, "submit-list": submit_list}[a.cmd](a)


if __name__ == "__main__":
    main()
