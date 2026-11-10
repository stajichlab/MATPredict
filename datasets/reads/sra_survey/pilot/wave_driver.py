#!/usr/bin/env python3
"""Run raw-read MAT typing for a WAVE of species: plan, test, gate, full run, report. Idempotent; every step can be repeated; state is files.

A wave directory (WAVE) holds:
  species.tsv                       input: species<TAB>slug (one row per species)
  plan/<slug>/                      panel/ (lineage panels L*/), runs_test.tsv, runs_full.tsv, qc_panel.json, status (planned | held:<reason>), thin (thin-evidence tokens, if any)
  out/                              ONE folder of per-run results for all species (<run>.<prefix>.tsv, <run>.run.json); a finished run is never redone
  test/, full/                      run lists, job ids and summaries of the two stages
  gates.tsv                         the decision per species and the numbers behind it
  WAVE_REPORT.md                    across species

Commands (all take --wave DIR and the config options below):
  plan-one SPECIES     build the lineage panels, check them (panel gate), choose test and full run lists for ONE species
  submit-plan          queue `plan-one` for every species of species.tsv as an array (panels are built on compute nodes)
  submit-test          merge the test runs of the planned species and queue them (split over the --parts partitions: stajichlab and preempt by default)
  gate                 summarize the test results and decide per species: PASS, PROVISIONAL, PROVISIONAL_HOMOTHALLIC, UNVALIDATED:<reasons>, NOT_TYPABLE:<cause> or HOLD:<reasons>   (gates below; causes: CAUSES)
  curate               apply the curator rulings (ZYGOSITY, RELEASED; see the code) to an existing gates table; a ZYGOSITY species is counted again from --stage's results
  incomplete           list the runs of --stage-name that did not finish (no run.json, `started` = killed or preempted, failed download/convert) into runs_incomplete.tsv for submit-list
  apply-causes         apply causes.tsv (curator ruling 2026-11-09: NOT_TYPABLE:homothallic / silent_cassettes, zygosity, HOLD:cause_unchecked) to the committed
                       wave folders (--pilot), writing cause_changes.tsv beside each changed gates_final.tsv
  regate               decide again from the numbers in an existing gates table (after the rules change); writes the previous decision beside the new one
  submit-full          queue the remaining runs of the species that are not held (PASS, PROVISIONAL, PROVISIONAL_HOMOTHALLIC and UNVALIDATED)
  report               write WAVE_REPORT.md from the summaries and the gates
  status               counts of species by state, and of runs done and failed

Gates. Panel (before any download): lineage panel L1 exists; at least 1 locus of each idiomorph went into the panel; core k-mers per idiomorph >= 300 and within a factor 10
of each other; shared flank k-mers >= 2000; at least 1 selectable run. Thin evidence (curator 2026-11-09): a species whose panel has only 1 locus of an idiomorph (single_locus,
for example a species with one assembled strain) or with fewer than 4 selectable runs (few_runs_<n>) is planned and typed, but its decision is capped at UNVALIDATED with
those reasons first (plan/<slug>/thin; ' thin=' in the gates table). Until 2026-11-09 both were plan-stage HOLDs (too_few_loci, too_few_runs_<n>). Test (after 20 runs): at least 90% of runs finished (runs not at NCBI are not counted); typed
(one idiomorph) >= 75% of finished runs and `both` <= 20%; median shared-flank depth >= 3; with >= 5 typed runs that have a `detect` truth: agreement >= 95% and at most one
disagreement; with fewer than 5 such runs the species is PROVISIONAL (typed >= 85% required). A species with at least half of its runs `both` and typed+both >= 85% is UNVALIDATED:cause_unchecked unless causes.tsv gives its cause (until 2026-11-09 it was PROVISIONAL_HOMOTHALLIC; only a cause with evidence now calls a species homothallic). Thresholds are in GATE below.

A species that has a usable panel and enough runs is not hidden for missing those limits (curator, 2026-11-02): it is UNVALIDATED when it types at least 50% of its finished runs and its only failures are the soft ones (typed
below 75% or 85%, `both` above 20%, shared depth below 3, or agreement short of the limit with fewer than 10 comparison runs). It gets the full run and is shown on the site as "unvalidated" with the reasons. HOLD is
left for what makes a call unusable: no panel, an idiomorph with no locus, or no runs (the plan status), fewer than half the finished runs typed, too many runs failing to download, or agreement short of the limit with 10 or more comparison runs.
"""
import argparse, collections, csv, json, os, re, statistics, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
P = "/bigdata/stajichlab/jstajich/projects/MATPredict"
PY = f"{P}/.pixi/envs/default/bin/python"
GATE = dict(min_loci=2, min_core_kmers=300, max_core_ratio=10, min_shared_kmers=2000, min_runs=4,
            min_done=0.90, min_typed=0.75, max_both=0.20, min_shared_depth=3.0, min_truth_runs=5, min_agree=0.95, max_disagree=1, provisional_typed=0.85, homothallic_both=0.5,
            unvalidated_typed=0.50, hold_truth_runs=10)
# min_loci and min_runs are thin-evidence floors (below them: UNVALIDATED with single_locus / few_runs_<n>); the hard floors are 1 locus per idiomorph and 1 run.
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
# Genera whose MAT loci fall outside the default span limits (4 to 60 kb): the Cryptococcus MAT locus is 40 to 200 kb (the a and alpha alleles span the whole region);
# the Rhodotorula P/R locus as reported is about 2 kb.
# Species complexes planned as ONE group (curator 2026-11-06): too few assemblies per species for a panel, enough together. The group name goes in species.tsv;
# its members (species-name prefixes) are what the panel builder and the run selector match.
GROUPS = {"Cryptococcus gattii complex": "Cryptococcus gattii|Cryptococcus deuterogattii|Cryptococcus bacillisporus",
          "Blumeria graminis complex": "Blumeria graminis|Blumeria hordei",                 # B. hordei was split from B. graminis (f. sp. hordei)
          "Sporothrix schenckii complex": "Sporothrix schenckii|Sporothrix brasiliensis"}
SPAN_GENERA = {"Cryptococcus": ("4000", "250000"), "Rhodotorula": ("1000", "60000"), "Rhizopus": ("1000", "60000"), "Mucor": ("1000", "60000"),
               "Ophidiomyces": ("4000", "20000")}      # Ophidiomyces: 49 of 50 MAT1-2 loci span ~11 kb; the 30 kb GCA_022829885.1 locus made a strain-specific MAT1-2 core (wave 1 untyped MAT1-2 runs, breadth ~0.4)
# Genera whose loci are reported as little more than their genes (no flank inside the reported span): cut this much genome on each side (panel_from_assemblies --pad).
PAD_GENERA = {"Rhodotorula": "6000", "Rhizopus": "6000", "Mucor": "6000"}


def plan_one(a):
    sp = a.species; sl = slug(sp); d = os.path.join(a.wave, "plan", sl); os.makedirs(d, exist_ok=True)
    I = a.inputs; prefix = sl.split("_")[0]; panel = os.path.join(d, "panel")
    log = open(os.path.join(d, "plan.log"), "w"); lenient = sp.split()[0] in LENIENT_GENERA
    def hold(reason):
        open(os.path.join(d, "status"), "w").write("held:" + reason + "\n"); print(sp, "HELD", reason); sys.exit(0)
    r = subprocess.run([PY, f"{HERE}/panel_from_assemblies.py", "--species", GROUPS.get(sp, sp), "--genomes", f"{I}/asco_genomes.tsv", "--loci", f"{I}/asco_loci.tsv", "--reports", f"{I}/reports_all.tar.zst",
                        "--genome-dir", f"{I}/genomes", "--out", panel, "--prefix", prefix, "--min-cluster-size", "1", "--max-reps", "6", "--per-lineage", "--jobs", "4"] + (["--min-span", "1000", "--classes", "mat_locus", "partial_locus", "idiomorph_gene_only"] if lenient else ["--min-span", SPAN_GENERA[sp.split()[0]][0], "--max-span", SPAN_GENERA[sp.split()[0]][1]] if sp.split()[0] in SPAN_GENERA else ["--min-span", "4000"]) + (["--borrow", sp.split()[0] + " "] if getattr(a, "borrow", False) or lenient else []) + (["--pad", PAD_GENERA[sp.split()[0]]] if sp.split()[0] in PAD_GENERA else []),
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
    thin = []      # thin evidence (curator 2026-11-09): the species goes on, capped at UNVALIDATED with the reason (plan/<slug>/thin; decide_species)
    if min(qc["loci_per_idiomorph"].values()) < 1: hold("too_few_loci")      # an idiomorph with no locus at all (not even borrowed): nothing to build its panel from
    if min(qc["loci_per_idiomorph"].values()) < GATE["min_loci"]: thin.append("single_locus")
    c = list(l1["core"].values())
    if min(c) < GATE["min_core_kmers"]: hold("core_kmers_low")
    if max(c) / max(1, min(c)) > GATE["max_core_ratio"]: hold("core_kmers_unbalanced")
    if l1["shared"] < GATE["min_shared_kmers"]: hold("shared_kmers_low")
    common = ["--species", GROUPS.get(sp, sp), "--group", sl, "--panel-dir", panel, "--panel-prefix", prefix, "--runs", f"{I}/fungi_wgs_illumina.tsv", "--meta", f"{I}/fungi_meta.tsv",
              "--names", *names] + (["--classes", "mat_locus", "partial_locus", "idiomorph_gene_only"] if lenient else []) + ["--studies", f"{I}/studies.tsv", "--assembly-summary", f"{I}/asm_summary_fungi.txt", "--genomes", f"{I}/asco_genomes.tsv", "--loci", f"{I}/asco_loci.tsv"]
    r = subprocess.run([PY, f"{HERE}/select_runs.py", *common, "--n-truth", "12", "--n-free", "8", "--out", f"{d}/runs_test.tsv"], capture_output=True, text=True); log.write(r.stdout + r.stderr)
    if r.returncode or not os.path.exists(f"{d}/runs_test.tsv"): hold("run_selection_failed")
    test = rd(f"{d}/runs_test.tsv")
    if not test: hold("too_few_runs_0")
    if len(test) < GATE["min_runs"]: thin.append(f"few_runs_{len(test)}")
    tf = os.path.join(d, "thin")
    if thin: open(tf, "w").write("+".join(thin) + "\n")
    elif os.path.exists(tf): os.remove(tf)
    r = subprocess.run([PY, f"{HERE}/select_runs.py", *common, "--n-truth", "25", "--n-free", "25", "--include", ",".join(x["run_accession"] for x in test), "--out", f"{d}/runs_full.tsv"], capture_output=True, text=True); log.write(r.stdout + r.stderr)
    open(os.path.join(d, "status"), "w").write("planned\n"); print(sp, "planned:", len(test), "test runs;", len(rd(f"{d}/runs_full.tsv")), "full runs")


# ------------------------------------------------------------------------------------------------ submit helpers
def sbatch(args, **kw):
    r = sh(["sbatch", "--parsable", *args]); 
    if r.returncode: sys.exit(f"sbatch failed: {r.stderr.strip()}")
    return r.stdout.strip().split(";")[0]


def parse_part(x):
    """'partition:account[:qos]' -> (partition, account, qos or ''); a bare name is its own account. The short partition needs its qos: short:stajichlab:short."""
    f = x.split(":")
    return (f[0], f[1] if len(f) > 1 else f[0], f[2] if len(f) > 2 else "")


INCOMPLETE = ("missing", "started", "failed_download", "failed_convert", "failed_unfinished")


def incomplete(a):
    """Runs of --stage-name STAGE that did not finish: no run.json (the task never started or died at once), state `started` (killed mid-run: preempted, node
    failure, time limit) or a failed download or conversion. Writes WAVE/STAGE/runs_incomplete.tsv (the stage's own run-list rows, ready for submit-list) and prints
    the count per state. `failed_type` and `unavailable_at_ncbi` are not rerun: the first is a typing error to look at, the second is not at NCBI."""
    d = os.path.join(a.wave, a.stage_name); out = a.out_dir or os.path.join(a.wave, "out")
    rows = rd(os.path.join(d, "runs_all.tsv")); bad = []; c = collections.Counter()
    for r in rows:
        j = os.path.join(out, r["run_accession"] + ".run.json")
        try:
            st = json.load(open(j)).get("state", "started") if os.path.exists(j) else "missing"
        except ValueError:
            st = "failed_unfinished"
        c[st] += 1
        if st in INCOMPLETE: bad.append(r)
    if bad: wr(os.path.join(d, "runs_incomplete.tsv"), bad, list(rows[0].keys()))
    print(f"{a.stage_name}: {dict(c)}; {len(bad)} to rerun" + (f" -> {d}/runs_incomplete.tsv" if bad else ""))


def submit_stage(a, stage, rows):
    """rows: run-list rows of ALL species for this stage; split over the two partitions; arrays of at most 2000 tasks; retry pass; summary job."""
    d = os.path.join(a.wave, stage); os.makedirs(d, exist_ok=True); os.makedirs(os.path.join(d, "logs"), exist_ok=True)
    out = a.out_dir or os.path.join(a.wave, "out"); os.makedirs(out, exist_ok=True)
    cols = list(rows[0].keys()); allp = os.path.join(d, "runs_all.tsv"); wr(allp, rows, cols)
    parts = [parse_part(x) for x in getattr(a, "parts", "stajichlab:stajichlab,preempt:preempt").split(",")]      # partition:account[:qos]; the runs are dealt over them in turn
    rpart = parse_part(getattr(a, "retry_part", None) or "stajichlab:stajichlab")      # the retry pass never runs on a preemptible partition: a run killed twice would be lost silently
    no_retry = rpart[0] == "none"                     # --retry-part none: no retry arrays (halves the queued tasks); rerun what `incomplete` lists afterwards
    spart = ("stajichlab", "stajichlab", "") if no_retry else rpart
    ids = []; jobs = []
    for pi, (part, acct, qos) in enumerate(parts):
        mine = rows[pi::len(parts)]
        for ci in range(0, len(mine), 2000):
            chunk = mine[ci:ci + 2000]; f = os.path.join(d, f"runs_{part}_{ci // 2000}.tsv"); wr(f, chunk, cols)
            for attempt in (("run",) if no_retry else ("run", "retry")):
                dep = ["--dependency", f"afterany:{ids[-1]}"] if attempt == "retry" else []
                p_, a_, q_ = (rpart if attempt == "retry" else (part, acct, qos))
                jid = sbatch(["-p", p_, "-A", a_, *(["--qos", q_] if q_ else []), "-c", str(a.cpus), "--mem", a.mem, "-J", f"wave-{stage}-{part}-{attempt}", "--array", f"1-{len(chunk)}%{a.conc}", *dep, "-o", f"{d}/logs/%x_%A_%a.out",
                              "--export", f"ALL,WORKER_DIR={HERE},TARGET_SPOTS=2000000,WINDOWS=8" + "".join("," + e for e in a.env), f"{HERE}/stream_type_sra.slurm", f, "none", "none", out]) if not a.dry_run else f"dry{len(ids)}"
                ids.append(jid); jobs.append((stage, p_, attempt, jid, len(chunk)))
    last = ids[-1]
    sm = sbatch(["-p", spart[0], "-A", spart[1], *(["--qos", spart[2]] if spart[2] else []), "-c", "1", "--mem", "4gb", "-t", "0:30:00", "-J", f"wave-{stage}-summary", "--dependency", "afterany:" + ":".join(ids), "-o", f"{d}/logs/%x_%j.out",
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
    jid = sbatch(["-p", "stajichlab", "-A", "stajichlab", "-c", "4", "--mem", "16gb", "-t", "2:00:00", "-J", "wave-plan", "--array", f"1-{len(sp)}%{a.conc}", "-o", f"{a.wave}/plan_logs/%x_%A_%a.out", "--wrap", cmd]) if not a.dry_run else "dryplan"
    print("plan array", jid, "for", len(sp), "species")


def submit_test(a):
    rows = []
    for sl in planned(a): rows += rd(os.path.join(a.wave, "plan", sl, "runs_test.tsv"))
    if not rows: sys.exit("no planned species")
    submit_stage(a, "test", rows)


# ------------------------------------------------------------------------------------------------ gate
# Curator rulings (2026-11-06) that the numbers alone cannot make.
# ZYGOSITY: a heterozygous diploid whose usual state is a/alpha (Candida albicans: most isolates are MTL-heterozygous). `both` is then a call (MTL a/alpha), so typed = one idiomorph
#   or both, the `both` limit does not apply, and an assembly's single idiomorph agrees with `both` (an assembly of a heterozygote collapses to one haplotype). Passing words carry _ZYGOSITY.
# RELEASED: species held only because the reads do not separate the active MAT locus from the silent HML/HMR cassettes (Saccharomyces; Ogataea's two flip-flop loci). The calls are not
#   wrong, the panel cannot tell the loci apart; they get the full run, flagged UNVALIDATED:released_<why>, so they can be examined (not when the runs failed to download).
ZYGOSITY = {"Candida albicans", "Pichia kudriavzevii"}      # P. kudriavzevii: most isolates are diploid a/alpha (curator 2026-11-06; wave 3 test 8 of 10 both)
RELEASED = {"Saccharomyces paradoxus": "cassettes", "Saccharomyces eubayanus": "cassettes", "Saccharomyces pastorianus": "cassettes", "Saccharomyces mikatae": "cassettes",
            "Ogataea haglerorum": "cassettes"}
# CAUSES (curator ruling 2026-11-09, #195/#196; design docs/notes/2026-11-09_junction-panel-design.md): why presence typing cannot give a species' mating type, from causes.tsv
# (species, cause, build_cause, basis; the evidence is cause_evidence.py's table). The cause overrides the numbers:
#   homothallic        both MAT genes in every genome (one locus or unlinked loci; the layout is in basis) -> NOT_TYPABLE:homothallic (no allele frequencies on the site)
#   silent_cassettes   expressed locus plus silent copies (HML/HMR, flip-flop, mat2/3) -> NOT_TYPABLE:silent_cassettes
#   zygosity           heterozygous diploid: joins ZYGOSITY (both = MTL a/alpha)
#   cause_unchecked    `both` in most runs for a reason not yet checked                -> UNVALIDATED:cause_unchecked+<soft reasons> (shown with counts and the
#                      warning; HOLD only for many failed downloads). A species with `both` in at least half its runs and NO cause gets the same (decide();
#                      it was PROVISIONAL_HOMOTHALLIC until 2026-11-09): only a causes.tsv row with evidence calls a species homothallic.
# A species with only a build_cause keeps its plan-stage HOLD (a technical reason). NOT_TYPABLE species are not queued again (submit-full skips them).
NOT_TYPABLE_CAUSES = ("homothallic", "silent_cassettes")
CAUSES_FILE = os.path.join(HERE, "causes.tsv")


def load_causes(path=CAUSES_FILE):
    """{species: cause} for the rows of causes.tsv with a cause (empty when the file is missing)."""
    return {r["species"]: r["cause"] for r in rd(path) if r.get("cause")} if os.path.exists(path) else {}


CAUSES = load_causes()
ZYGOSITY |= {s for s, c in CAUSES.items() if c == "zygosity"}


def cause_decision(species):
    """The decision a cause forces whatever the numbers (NOT_TYPABLE:<cause>), or None (no cause, zygosity, or cause_unchecked, which is decided from the numbers)."""
    c = CAUSES.get(species)
    return f"NOT_TYPABLE:{c}" if c in NOT_TYPABLE_CAUSES else None


THIN_RE = re.compile(r"\bthin=(\S+)")


def thin_of(reasons):
    """The thin-evidence tokens recorded in a gates-table reasons text (' thin=single_locus+few_runs_2'), as a tuple."""
    m = THIN_RE.search(reasons or "")
    return tuple(m[1].split("+")) if m else ()


def plan_thin(wave, sl):
    """The thin-evidence tokens plan-one wrote to plan/<slug>/thin (single_locus, few_runs_<n>), as a tuple."""
    f = os.path.join(wave, "plan", sl, "thin")
    return tuple(t for t in open(f).read().strip().split("+") if t) if os.path.exists(f) else ()


def decide_species(species, n_counted, n_done, n_typed, n_both, n_truth, n_agree, msd, thin=()):
    """decide() with the curator rulings for SPECIES applied. For a ZYGOSITY species n_typed must already include the `both` runs and n_agree count `both` as agreeing (numbers()).
    thin: thin-evidence tokens from the plan (single_locus: an idiomorph's panel comes from one assembly locus; few_runs_<n>: fewer than GATE min_runs runs exist). Curator
    2026-11-09: such a species is typed and shown, but never above UNVALIDATED, with the tokens first among the reasons; a HOLD stays a HOLD (the tokens are added)."""
    dec, why = _decide_species(species, n_counted, n_done, n_typed, n_both, n_truth, n_agree, msd)
    if not thin or dec.startswith("NOT_TYPABLE"): return dec, why
    why = list(thin) + [w for w in why if w not in thin]
    return ("HOLD:" if dec.startswith("HOLD") else "UNVALIDATED:") + "+".join(why), why


def _decide_species(species, n_counted, n_done, n_typed, n_both, n_truth, n_agree, msd):
    forced = cause_decision(species)
    if forced: return forced, [forced.split(":", 1)[1]]
    if species in ZYGOSITY:
        dec, why = decide(n_counted, n_done, n_typed, 0, n_truth, n_agree, msd, zygosity=True)
    else:
        dec, why = decide(n_counted, n_done, n_typed, n_both, n_truth, n_agree, msd)
    if species in RELEASED and dec.startswith("HOLD") and "many_failed" not in why:
        why = [f"released_{RELEASED[species]}"] + why
        dec = "UNVALIDATED:" + "+".join(why)
    if CAUSES.get(species) == "cause_unchecked" and "many_failed" not in why:      # curator 2026-11-09: kept visible with the warning
        why = ["cause_unchecked"] + [w for w in why if w != "cause_unchecked"]
        dec = "UNVALIDATED:" + "+".join(why)
    return dec, why


def numbers(rs, zygosity=False):
    """(counted, done, typed, both, truth_typed, agree, median shared depth) from the result rows of one species. zygosity: `both` is a call and agrees with either single truth."""
    counted = [x for x in rs if x["state"] != "unavailable_at_ncbi"]
    done = [x for x in counted if x["state"] == "done"]
    single = [x for x in done if x["call"] not in ("both", "none", "low_depth", "")]
    both = [x for x in done if x["call"] == "both"]
    typed = single + both if zygosity else single
    tr = [x for x in typed if x["truth"] and "," not in x["truth"]]
    agree = sum(1 for x in tr if x["call"] == x["truth"] or (zygosity and x["call"] == "both"))
    sd = [float(x["shared_depth"]) for x in done if x["shared_depth"] not in ("", None)]
    return len(counted), len(done), len(typed), (0 if zygosity else len(both)), len(tr), agree, (statistics.median(sd) if sd else 0)


def decide(n_counted, n_done, n_typed, n_both, n_truth, n_agree, msd, zygosity=False):
    """(decision, reasons) for one tested species from the numbers of its test runs: runs counted (not "unavailable at NCBI"), finished, typed
    (one idiomorph), both, typed runs with a `detect` truth, those that agree, and the median shared-flank depth. A reason is "hard" when the
    calls cannot be used at all (HOLD) and "soft" when they can be used with the reason beside them (UNVALIDATED)."""
    hard, soft = [], []
    tf = n_typed / max(1, n_done); bf = n_both / max(1, n_done)
    if not n_counted or n_done / n_counted < GATE["min_done"]: hard.append("many_failed")
    homo = not zygosity and bf >= GATE["homothallic_both"] and (n_typed + n_both) / max(1, n_done) >= GATE["provisional_typed"]    # a homothallic species: `both` is the expected answer
    if not homo:
        if tf < GATE["min_typed"]: (hard if tf < GATE["unvalidated_typed"] else soft).append(f"typed_{tf:.2f}")
        if bf > GATE["max_both"]: soft.append(f"both_{bf:.2f}")
    if msd < GATE["min_shared_depth"]: soft.append(f"shared_depth_{msd:.1f}")
    prov = n_truth < GATE["min_truth_runs"] or homo
    if not prov:
        if n_agree / n_truth < GATE["min_agree"] or n_truth - n_agree > GATE["max_disagree"]:
            (hard if n_truth >= GATE["hold_truth_runs"] else soft).append(f"agree_{n_agree}of{n_truth}")
    elif not homo and tf < GATE["provisional_typed"]:
        (hard if tf < GATE["unvalidated_typed"] else soft).append(f"provisional_typed_{tf:.2f}")
    if zygosity and (hard or soft): soft = ["zygosity"] + soft
    if homo and not hard:                            # mostly `both` with no cause on file: no claim of homothallism (curator 2026-11-09, #195)
        why = ["cause_unchecked", f"both_{bf:.2f}"] + soft
        return "UNVALIDATED:" + "+".join(why), why
    if hard: return "HOLD:" + "+".join(hard + soft), hard + soft
    if soft: return "UNVALIDATED:" + "+".join(soft), soft
    return ("PROVISIONAL" if prov else "PASS") + ("_ZYGOSITY" if zygosity else ""), []


def gate_row(sl, species, rs, thin=()):
    """One gates-table row from the result rows of a species. The reasons text holds the numbers regate() reads; a ZYGOSITY species' numbers are already counted that way (mode=zygosity);
    thin-evidence tokens from the plan are kept as ' thin=a+b'."""
    zyg = species in ZYGOSITY
    n_counted, n_done, n_typed, n_both, n_truth, n_agree, msd = numbers(rs, zyg)
    dec, _ = decide_species(species, n_counted, n_done, n_typed, n_both, n_truth, n_agree, msd, thin)
    return dict(slug=sl, species=species, decision=dec, reasons=f"runs={len(rs)} done={n_done} typed={n_typed} both={n_both} none/low={n_done - n_typed - n_both} truth_typed={n_truth} agree={n_agree} median_shared_depth={msd:.1f}" + (" mode=zygosity" if zyg else "") + (" thin=" + "+".join(thin) if thin else ""))


def gate(a):
    res = rd(os.path.join(a.wave, a.stage, "summary", "results_all.tsv"))
    by = collections.defaultdict(list)
    for r in res: by[r["group"]].append(r)
    rows = []
    for r0 in rd(os.path.join(a.wave, "species.tsv")):
        sl = r0["slug"]; sf = os.path.join(a.wave, "plan", sl, "status"); st = open(sf).read().strip() if os.path.exists(sf) else "not_planned"
        if st != "planned":
            rows.append(dict(slug=sl, species=r0["species"], decision=cause_decision(r0["species"]) or "HOLD:" + st.replace("held:", ""), reasons=st)); continue
        rs = by.get(sl, [])
        if not rs: rows.append(dict(slug=sl, species=r0["species"], decision="PENDING", reasons="no test results yet")); continue
        rows.append(gate_row(sl, r0["species"], rs, plan_thin(a.wave, sl)))
    wr(os.path.join(a.wave, "gates.tsv" if a.stage == "test" else f"gates_{a.stage}.tsv"), rows, ["slug", "species", "decision", "reasons"])
    print(collections.Counter(r["decision"].split(":")[0] for r in rows))


REGATE = re.compile(r"runs=(\d+) done=(\d+) typed=(\d+) both=(\d+) none/low=(\d+) truth_typed=(\d+) agree=(\d+) median_shared_depth=([\d.]+)")


def regate(a):
    """Decide again from the numbers a gates table already holds (its `reasons` column), for example after the rules changed. A species
    that was not tested (a plan status such as held:panel_build_failed) keeps its decision; one held for many failed downloads stays held.
    Only a species whose decision changes is rewritten (a species held before and now keeps its row). Writes --out-file (default: the input file,
    in place) and, beside it, regate_changes.tsv with the previous and the new decision of each species that changed."""
    src = os.path.join(a.wave, a.in_file); dst = os.path.join(a.wave, a.out_file or a.in_file)
    rows = rd(src); changes = []
    for r in rows:
        m = REGATE.search(r["reasons"])
        if not m:
            forced = cause_decision(r["species"])             # a species held at the panel stage still takes its cause
            if forced and forced != r["decision"]:
                changes.append({"slug": r["slug"], "species": r["species"], "previous": r["decision"], "decision": forced, "reasons": r["reasons"]}); r["decision"] = forced
            continue
        if r["species"] in ZYGOSITY and "mode=zygosity" not in r["reasons"]:
            continue                                                                   # counted the ordinary way: only `curate` (from the results) can decide it again
        done, typed, both, truth, agree, msd = int(m[2]), int(m[3]), int(m[4]), int(m[6]), int(m[7]), float(m[8])
        counted = done if "many_failed" not in r["decision"] else done * 2 + 1       # the counted runs are not in the table; this keeps a species held for failed downloads held
        thin = thin_of(r["reasons"])
        new, _ = decide_species(r["species"], counted, done, typed, both, truth, agree, msd, thin)
        if new.split(":")[0] != r["decision"].split(":")[0] or ("released_" in new) != ("released_" in r["decision"]) or (new != r["decision"] and ("NOT_TYPABLE" in new or "cause_unchecked" in new or thin)):
            changes.append({"slug": r["slug"], "species": r["species"], "previous": r["decision"], "decision": new, "reasons": r["reasons"]})
            r["decision"] = new
    wr(dst, rows, ["slug", "species", "decision", "reasons"])
    wr(os.path.join(os.path.dirname(dst), "regate_changes.tsv"), changes, ["slug", "species", "previous", "decision", "reasons"])
    kinds = collections.Counter((c["previous"].split(":")[0], c["decision"].split(":")[0]) for c in changes)
    print(f"{len(rows)} species -> {dst}; decisions: {dict(collections.Counter(r['decision'].split(':')[0] for r in rows))}; changed: {dict(kinds)}")


def curate(a):
    """Apply the curator rulings (ZYGOSITY, RELEASED) to an existing gates table (--in-file, in place or to --out-file). A ZYGOSITY species is counted again from the results of
    --stage (it needs the calls, not only the numbers); a RELEASED species is decided again from its numbers. Other species are untouched. Changes go to curate_changes.tsv."""
    src = os.path.join(a.wave, a.in_file); dst = os.path.join(a.wave, a.out_file or a.in_file)
    rows = rd(src); changes = []; res = None
    for i, r in enumerate(rows):
        if r["species"] not in ZYGOSITY and r["species"] not in RELEASED and r["species"] not in CAUSES:
            continue
        if cause_decision(r["species"]) and not REGATE.search(r["reasons"]):
            new = dict(r, decision=cause_decision(r["species"]))                      # held at the panel stage: the cause alone decides
        elif r["species"] in ZYGOSITY:
            if res is None:
                res = collections.defaultdict(list)
                for x in rd(os.path.join(a.wave, a.stage, "summary", "results_all.tsv")): res[x["group"]].append(x)
            if not res.get(r["slug"]):
                continue
            new = gate_row(r["slug"], r["species"], res[r["slug"]], thin_of(r["reasons"]))
        else:
            m = REGATE.search(r["reasons"])
            if not m:
                continue
            done = int(m[2]); counted = done if "many_failed" not in r["decision"] else done * 2 + 1
            new = dict(r, decision=decide_species(r["species"], counted, done, int(m[3]), int(m[4]), int(m[6]), int(m[7]), float(m[8]), thin_of(r["reasons"]))[0])
        if new != r:
            changes.append({"slug": r["slug"], "species": r["species"], "previous": r["decision"], "decision": new["decision"], "reasons": new["reasons"]}); rows[i] = new
    wr(dst, rows, ["slug", "species", "decision", "reasons"])
    wr(os.path.join(os.path.dirname(dst), "curate_changes.tsv"), changes, ["slug", "species", "previous", "decision", "reasons"])
    for c in changes: print(f"{c['species']}: {c['previous']} -> {c['decision']}")
    print(len(changes), "species changed ->", dst)


def apply_causes(a):
    """Apply causes.tsv (and the ZYGOSITY counting it adds) to the committed wave folders under --pilot: in each wave*/gates_final.tsv, a species with a cause is
    decided again, from its rows in results_all_final.tsv when it was typed (gate_row: a zygosity species is counted that way) and from the cause alone when the panel
    stage held it. Rewrites gates_final.tsv in place and writes cause_changes.tsv beside it (slug, species, previous, decision, reasons) for every changed row."""
    total = 0
    for w in sorted(os.listdir(a.pilot)):
        wd = os.path.join(a.pilot, w); gf = os.path.join(wd, "gates_final.tsv"); rf = os.path.join(wd, "results_all_final.tsv")
        if not (re.fullmatch(r"wave\d+[a-z_]*", w) and os.path.exists(gf)): continue
        rows = rd(gf); res = collections.defaultdict(list); changes = []
        if os.path.exists(rf):
            for x in rd(rf): res[x["group"]].append(x)
        for i, r in enumerate(rows):
            if r["species"] not in CAUSES: continue
            if r["reasons"].startswith("runs=") and res.get(r["slug"]): new = gate_row(r["slug"], r["species"], res[r["slug"]], thin_of(r["reasons"]))
            else: new = dict(r, decision=cause_decision(r["species"]) or r["decision"])
            if new["decision"] != r["decision"] or new["reasons"] != r["reasons"]:
                changes.append({"slug": r["slug"], "species": r["species"], "previous": r["decision"], "decision": new["decision"], "reasons": new["reasons"]}); rows[i] = new
        if changes:
            wr(gf, rows, ["slug", "species", "decision", "reasons"])
            wr(os.path.join(wd, "cause_changes.tsv"), changes, ["slug", "species", "previous", "decision", "reasons"])
        for c in changes: print(f"{w}\t{c['species']}\t{c['previous']} -> {c['decision']}")
        total += len(changes)
    print(total, "species rows changed")


def submit_full(a):
    """Queue the runs of the species that are not held. --only WORD limits it to the species whose decision is WORD (for example UNVALIDATED), --stage-name
    names the stage folder (default full; a second pass must not overwrite the first), --gates the gates table, and --out-dir a results folder of its own (a run whose result is
    already in the default folder is skipped by the worker, so a pass that adds a screen needs its own)."""
    g = {r["slug"]: r["decision"] for r in rd(os.path.join(a.wave, a.gates))}
    rows = []
    for sl, dec in g.items():
        word = dec.split(":")[0]
        if word in ("PASS", "PROVISIONAL", "PROVISIONAL_HOMOTHALLIC", "PASS_ZYGOSITY", "PROVISIONAL_ZYGOSITY", "UNVALIDATED") and (not a.only or word == a.only):
            rows += rd(os.path.join(a.wave, "plan", sl, "runs_full.tsv"))[:a.cap]     # the test runs come first; at most --cap runs per species
    if not rows: sys.exit("no species to run" + (f" with the decision {a.only}" if a.only else " (none passed the gate)"))
    submit_stage(a, a.stage_name if a.stage_name != "list" else "full", rows)


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
    ap.add_argument("cmd", choices=["plan-one", "submit-plan", "submit-test", "gate", "regate", "curate", "incomplete", "submit-full", "report", "status", "submit-list", "apply-causes"])
    ap.add_argument("--wave", help="the wave directory (every command but apply-causes)"); ap.add_argument("--pilot", default=HERE, help="apply-causes: the folder of committed wave folders (default: this script's folder)"); ap.add_argument("--inputs", default=f"{P}/results/sra_reads_typing/wave_inputs")
    ap.add_argument("--species"); ap.add_argument("--conc", type=int, default=60, help="tasks running at once per array (the curator lifted the CPU limit for this work, 2026-10-27)")
    ap.add_argument("--cpus", type=int, default=1, help="CPUs per run task: every step is single-threaded (measured 2026-10-27), and the exfab account is limited to 32 CPUs, so 1 CPU gives 32 tasks at once instead of 8")
    ap.add_argument("--mem", default="6gb", help="memory per run task (peak measured 4.9 GB)")
    ap.add_argument("--cap", type=int, default=50, help="most runs per species in the full stage")
    ap.add_argument("--stage", default="test", help="gate: which stage's summary to judge (a folder under the wave directory with summary/results_all.tsv) (the full stage has more runs per species)")
    ap.add_argument("--env", action="append", default=[], help="extra KEY=VALUE for the run tasks (for example READS_K=21); repeatable")
    ap.add_argument("--runs", help="submit-list: a run-list TSV with panel_dir and panel_prefix columns")
    ap.add_argument("--stage-name", default="list", help="submit-list and submit-full: stage folder name under the wave directory (submit-full defaults to full)")
    ap.add_argument("--out-dir", help="submit-list: results folder (default <wave>/out; the worker skips a run whose result is already there, so a different analysis of the same runs, for example sketches, needs its own)")
    ap.add_argument("--retry-held", action="store_true", help="submit-plan: only the species held for too few runs")
    ap.add_argument("--retry-status", help="submit-plan: only the species whose plan status starts with this (for example held:panel_build_failed)")
    ap.add_argument("--borrow", action="store_true", help="plan-one, submit-plan: an idiomorph with no locus in the species takes its loci from the other species of the genus")
    ap.add_argument("--gates", default="gates.tsv", help="submit-full: the gates table to read (in the wave directory; gates_final.tsv has the decisions after the retries)")
    ap.add_argument("--only", help="submit-full: only the species whose decision is this word (PASS, PROVISIONAL, PROVISIONAL_HOMOTHALLIC or UNVALIDATED)")
    ap.add_argument("--in-file", default="gates_final.tsv", help="regate: the gates table (in the wave directory) to decide again")
    ap.add_argument("--out-file", help="regate: where to write it (default: the input file)")
    ap.add_argument("--parts", default="stajichlab:stajichlab,preempt:preempt", help="partition:account pairs the run arrays are split over (the exfab partition is one node, gpu12, whose SRA downloads run at about 0.3 MB/s: not a default) (preempt:preempt is allowed: a preempted task is rerun by the retry pass, and the worker is idempotent)")
    ap.add_argument("--retry-part", default="stajichlab:stajichlab", help="partition:account[:qos] of the retry pass and the summary job; keep it NOT preemptible. 'none': no retry arrays (then run `incomplete` and resubmit)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if a.cmd != "apply-causes" and not a.wave: ap.error("--wave is required")
    {"apply-causes": apply_causes, "plan-one": plan_one, "submit-plan": submit_plan, "submit-test": submit_test, "gate": gate, "regate": regate, "curate": curate, "incomplete": incomplete, "submit-full": submit_full, "report": report, "status": status, "submit-list": submit_list}[a.cmd](a)


if __name__ == "__main__":
    main()
