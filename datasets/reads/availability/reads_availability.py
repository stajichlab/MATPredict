#!/usr/bin/env python3
"""Which public and local read sets could still be MAT-typed, per species: PopSeq's catalog joined to what MATPredict has typed.

    reads_availability.py --popseq /bigdata/stajichlab/shared/projects/Fungi/PopSeq \
        --local-inventory INVENTORY_OUT --db ../../../db --out snapshots/ [--previous snapshots/2026-10-26]

Sources (nothing is downloaded here):
- PopSeq (the lab's shared catalog) says what EXISTS: the newest kingdom run table `cache/fungi_runs_ena_<date>.tsv.gz`
  (ENA read_run, Illumina WGS GENOMIC under tax_tree(4751)), `catalog/taxonomy_species_map.csv.gz` (any taxid -> its species taxid
  and current NCBI name) and `catalog/genome_sizes.csv` (per-species genome size, for a coverage estimate).
- MATPredict says what has been TYPED: the committed reads result tables (`--typed`, `run` + `call` + `state`), the wave gates
  (`--gates`) and the Coccidioides local-reads table (`--typed-local`, keyed by strain).
- PopSeq's local inventory (`popseq_inventory.py --out-root DIR`, spec 005) says what the lab holds on disk under
  Population_Genomics; `local_subprojects.tsv` maps each subproject to a species.

Eligibility repeats the survey's rule (survey_ena.py): a study of the species that classify_study() does not flag
(mutant, evolution, cross, not_wgs) and a run of at least 300 Mb, paired or single-end. One isolate is one BioSample.
Writes OUT/<run-table date>/: species_availability.tsv, untyped_runs.tsv.gz, local_samples.tsv, SUMMARY.md, snapshot.json,
and CHANGES.md when --previous is given. Deterministic: same inputs, same bytes (snapshot.json records the inputs and their sha256).
Standard library only, except --db (needs PyYAML, as survey_ena.py).
"""
import argparse, collections, csv, glob, gzip, hashlib, json, math, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sra_survey"))
from survey_ena import classify_study, species as parse_species, NOT_SPECIES, MIN_BASES, db_tiers, tier  # noqa: E402

PILOT = os.path.join(HERE, "..", "sra_survey", "pilot")
WAVES = sorted(d for d in glob.glob(os.path.join(PILOT, "wave*")) if re.fullmatch(r"wave\d+[a-z_]*", os.path.basename(d)) and os.path.exists(os.path.join(d, "gates_final.tsv")))      # every committed wave (wave1, wave2, wave3_*, wave4*), as the reads site
DEFAULT_TYPED = [os.path.join(w, "results_all_final.tsv") for w in WAVES] + [os.path.join(PILOT, "pilot4", "results_all.tsv")]
DEFAULT_GATES = [os.path.join(w, "gates_final.tsv") for w in WAVES]
sys.path.insert(0, PILOT)
from wave_driver import GROUPS  # noqa: E402  (species complexes typed as one group; standard library only)
DEFAULT_TYPED_LOCAL =[os.path.join(HERE, "..", "coccidioides", "results", "cocci2_reads_vs_detect.tsv")]
NOT_A_CALL = {"", "both", "none", "low_depth", "untyped"}
PASSING = ("PASS", "PROVISIONAL")                      # PROVISIONAL_HOMOTHALLIC is handled separately
PANEL_HOLD = re.compile(r"panel_build_failed|too_few_loci|panel_L1_error|core_kmers_low|core_kmers_unbalanced|shared_kmers_low|panel_no_unique_kmers|run_selection_failed")
TARGET_BASES, MIN_SEL, MAX_SEL = 3e9, 0.6e9, 12e9     # select_runs.py: nearest 3 Gb within 0.6-12 Gb


def opener(p):
    return gzip.open(p, "rt", newline="") if p.endswith(".gz") else open(p, newline="")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def norm_strain(s):
    """Strain names for a weak cross-source match: lower case, punctuation and spaces dropped."""
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def wilson_halfwidth(k, n):
    if n == 0:
        return ""
    z, p = 1.96, k / n
    return f"{z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n):.3f}"


def newest_run_table(root):
    c = sorted(glob.glob(os.path.join(root, "cache", "fungi_runs_ena_*.tsv.gz")))
    if not c:
        sys.exit(f"no cache/fungi_runs_ena_*.tsv.gz under {root}; pass --runs")
    return c[-1]


def load_taxmap(path):
    """taxid -> (species_taxid, species_name) for taxids at or below species; name -> species name for the same rows."""
    by_id, by_name = {}, {}
    with opener(path) as fh:
        for r in csv.DictReader(fh):
            if r["resolution"].startswith("above_species") or not r["species_taxid"]:
                continue
            by_id[r["taxid"]] = (r["species_taxid"], r["species_name"])
            by_name.setdefault(r["name"], (r["species_taxid"], r["species_name"]))
    return by_id, by_name


def resolve(name, by_name):
    """A name from a MATPredict table (deposited spelling) -> (species_taxid, species_name); falls back to the name itself."""
    if name in by_name:
        return by_name[name]
    return ("", name)


def load_runs(path, by_id, skip):
    """Per species: runs (dicts) of the kingdom table that name a species, are Illumina WGS and at least MIN_BASES.
    Returns ({species_taxid: {name, aliases, runs: [..]}}, run -> species_taxid for every run, run -> biosample, counts)."""
    sp = {}
    run_sp, run_bs, n = {}, {}, collections.Counter()
    with opener(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n["runs"] += 1
            res = by_id.get(r["tax_id"])
            if res is None:
                name = parse_species(r["scientific_name"])
                if not name or NOT_SPECIES.search(r["scientific_name"]):
                    n["no_species"] += 1
                    continue
                res = ("name:" + name, name)
            key, name = res
            run_sp[r["run_accession"]] = key
            run_bs[r["run_accession"]] = r["sample_accession"].split(";")[0]
            if name in skip:
                n["skipped_species"] += 1
                continue
            try:
                bases = int(r["base_count"] or 0)
            except ValueError:
                bases = 0
            if bases < MIN_BASES:
                n["below_min_bases"] += 1
                continue
            d = sp.setdefault(key, {"name": name, "aliases": set(), "runs": []})
            d["aliases"].add(" ".join(r["scientific_name"].split()[:2]))
            r["bases"] = bases
            d["runs"].append(r)
            n["kept"] += 1
    return sp, run_sp, run_bs, n


def eligible_isolates(runs):
    """Apply the survey's study flags; one isolate per BioSample (its run nearest 3 Gb inside 0.6-12 Gb, else the largest).
    Returns ({biosample: run}, {flag: isolates excluded}, {study: isolates kept})."""
    by_study = collections.defaultdict(list)
    for r in runs:
        by_study[r["study_accession"]].append(r)
    keep, excluded, per_study = {}, collections.Counter(), collections.Counter()
    for st in sorted(by_study):
        rs = by_study[st]
        texts = {}
        for r in rs:
            k = (r.get("strain") or r.get("isolate") or r["sample_accession"] or r["run_accession"]).strip().lower()
            texts[k] = " ".join(r.get(f, "") for f in ("sample_title", "strain", "isolate"))
        title = rs[0].get("study_title", "")
        flag = classify_study(title, texts)
        bss = {r["sample_accession"].split(";")[0] or r["run_accession"] for r in rs}
        if flag:
            excluded[flag] += len(bss)
            continue
        for r in rs:
            bs = r["sample_accession"].split(";")[0] or r["run_accession"]
            prev = keep.get(bs)
            if prev is None or rank(r) < rank(prev):
                keep[bs] = r
    for bs, r in keep.items():
        per_study[r["study_accession"]] += 1
    return keep, excluded, per_study


def rank(r):
    b = r["bases"]
    inside = MIN_SEL <= b <= MAX_SEL
    return (0 if inside else 1, abs(b - TARGET_BASES) if inside else -b, r["run_accession"])


def read_tsv(p):
    with opener(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--popseq", help="PopSeq root; supplies --runs, --taxmap and --genome-sizes when they are not given")
    ap.add_argument("--runs", help="PopSeq kingdom run table (default: the newest ROOT/cache/fungi_runs_ena_*.tsv.gz)")
    ap.add_argument("--taxmap", help="default ROOT/catalog/taxonomy_species_map.csv.gz")
    ap.add_argument("--genome-sizes", help="default ROOT/catalog/genome_sizes.csv")
    ap.add_argument("--popseq-commit", default="", help="PopSeq git commit the tables came from (read with git when ROOT is a checkout)")
    ap.add_argument("--typed", nargs="*", default=DEFAULT_TYPED, help="MATPredict reads result tables (run, call, state)")
    ap.add_argument("--gates", nargs="*", default=DEFAULT_GATES)
    ap.add_argument("--typed-local", nargs="*", default=DEFAULT_TYPED_LOCAL, help="tables of reads typed from local files (strain, species, reads_call)")
    ap.add_argument("--local-inventory", help="out-root of PopSeq popseq_inventory.py (holds panels/<subproject>/metadata.csv)")
    ap.add_argument("--local-map", default=os.path.join(HERE, "local_subprojects.tsv"))
    ap.add_argument("--db", help="MATPredict db folder: reference tier per species (needs PyYAML)")
    ap.add_argument("--skip-species", default="Saccharomyces cerevisiae")
    ap.add_argument("--min-isolates", type=int, default=4)
    ap.add_argument("--previous", help="an earlier snapshot folder; writes CHANGES.md against it")
    ap.add_argument("--out", required=True, help="snapshot root; the snapshot goes in OUT/<run table date>/")
    a = ap.parse_args()
    root = a.popseq
    if not (a.runs and a.taxmap) and not root:
        sys.exit("give --popseq, or --runs and --taxmap")
    a.runs = a.runs or newest_run_table(root)
    a.taxmap = a.taxmap or os.path.join(root, "catalog", "taxonomy_species_map.csv.gz")
    if not a.genome_sizes and root and os.path.exists(os.path.join(root, "catalog", "genome_sizes.csv")):
        a.genome_sizes = os.path.join(root, "catalog", "genome_sizes.csv")
    if not a.popseq_commit and root and os.path.isdir(os.path.join(root, ".git")):
        a.popseq_commit = subprocess.run(["git", "-C", root, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    m = re.search(r"(\d{4}-\d{2}-\d{2})", os.path.basename(a.runs))
    date = m.group(1) if m else "undated"
    out = os.path.join(a.out, date)
    os.makedirs(out, exist_ok=True)

    skip = {x.strip() for x in a.skip_species.split(",") if x.strip()}
    by_id, by_name = load_taxmap(a.taxmap)
    sp, run_sp, run_bs, counts = load_runs(a.runs, by_id, skip)
    gsize = {}
    if a.genome_sizes:
        for r in csv.DictReader(open(a.genome_sizes, newline="")):
            try:
                gsize[r["species_taxid"]] = float(r["genome_mb"])
            except (KeyError, ValueError):
                pass

    # --- what MATPredict has typed -------------------------------------------------------------
    typed = collections.defaultdict(dict)        # species key -> biosample (or run) -> call
    typed_unplaced = collections.Counter()
    slug_species = {}
    for p in a.typed:
        for r in read_tsv(p):
            if r.get("state") != "done":
                continue
            run = r["run"]
            key = run_sp.get(run)
            if key is None:
                typed_unplaced[r.get("group", "")] += 1
                continue
            typed[key][run_bs.get(run) or run] = r.get("call", "")
            slug_species.setdefault(r.get("group", ""), key)
    typed_local = collections.defaultdict(dict)  # species key -> normalised strain -> call
    for p in a.typed_local:
        for r in read_tsv(p):
            name = r["species"]
            key, _ = resolve(name, by_name)
            typed_local[key or "name:" + name][norm_strain(r["strain"])] = r.get("reads_call", "")
    gates = {}                                   # species key -> (decision, reasons, wave, slug in the wave)
    for p in a.gates:
        wave = os.path.basename(os.path.dirname(p))
        for r in read_tsv(p):
            members = GROUPS.get(r["species"])
            if members:          # a species complex typed as one group: every member species takes the group's decision and slug (its runs use the group panel)
                for m in members.split("|"):
                    mk = resolve(m, by_name)[0]
                    if mk: gates[mk] = (r["decision"], r["reasons"], wave, r["slug"])
                continue
            key = slug_species.get(r["slug"]) or resolve(r["species"], by_name)[0] or "name:" + r["species"]
            if key in gates and r["decision"].startswith("HOLD") and not gates[key][0].startswith("HOLD"):
                continue         # a member species' own held row (too few loci) must not undo its group's decision
            gates[key] = (r["decision"], r["reasons"], wave, r["slug"])

    # --- what the lab holds locally ------------------------------------------------------------
    local_rows, local_by_sp = [], collections.defaultdict(list)
    if a.local_inventory:
        lmap = {r["subproject"]: r for r in csv.DictReader(open(a.local_map, newline=""), delimiter="\t")}
        for meta in sorted(glob.glob(os.path.join(a.local_inventory, "panels", "*", "metadata.csv"))):
            sub = os.path.basename(os.path.dirname(meta))
            ent = lmap.get(sub) or {}
            name = ent.get("species", "")
            counted = name and ent.get("population", "yes") == "yes"     # a lab-lineage resequencing set is listed, not counted
            key = (resolve(name, by_name)[0] or "name:" + name) if name else ""
            for r in csv.DictReader(open(meta, newline="")):
                row = {"subproject": sub, "species": name, "species_taxid": key if not key.startswith("name:") else "",
                       "counted": "yes" if counted else "no", "sample": r["sample"], "provenance": r["provenance"],
                       "run_accession": r["run_accession"], "biosample": r["biosample"] or run_bs.get(r["run_accession"], ""), "n_files": r["n_files"]}
                local_rows.append(row)
                if counted:
                    local_by_sp[key].append(row)

    # --- per species ---------------------------------------------------------------------------
    spdb, gendb = db_tiers(a.db) if a.db else ({}, {})
    rows, untyped_rows = [], []
    keys = set(sp) | set(typed) | set(typed_local) | set(gates) | set(local_by_sp)
    for key in keys:
        d = sp.get(key, {"name": key[5:] if key.startswith("name:") else key, "aliases": set(), "runs": []})
        name = d["name"]
        if name in skip:
            continue
        elig, excl, per_study = eligible_isolates(d["runs"])
        public_names = {norm_strain(r.get("strain") or r.get("isolate")): bs for bs, r in elig.items() if (r.get("strain") or r.get("isolate"))}
        tp = typed.get(key, {})
        tl = typed_local.get(key, {})
        # local-typed strains that name a public isolate are counted as typed for that isolate (a weak, name-only match)
        tl_bs = {public_names[s] for s in tl if s in public_names}
        typed_bs = set(tp) | tl_bs
        untyped = [bs for bs in elig if bs not in typed_bs]
        calls = collections.Counter(tp.values()) + collections.Counter(tl.values())
        single = {c: n for c, n in calls.items() if c not in NOT_A_CALL}
        n_single = sum(single.values())
        attempted = sum(calls.values())
        minor = min(single.values()) if len(single) >= 2 else 0
        typed_rate = n_single / attempted if attempted else 0
        n_all = n_single + round(len(untyped) * typed_rate)
        gate = gates.get(key, ("", "", "", ""))
        decision = gate[0]
        loc = local_by_sp.get(key, [])
        seen = set()                         # one strain listed by two subprojects counts once
        loc_files = [r for r in loc if int(r["n_files"] or 0) > 0 and not (norm_strain(r["sample"]) in seen or seen.add(norm_strain(r["sample"])))]
        loc_public = [r for r in loc_files if r["run_accession"] or r["biosample"] in elig or norm_strain(r["sample"]) in public_names]
        loc_new = len(loc_files) - len(loc_public)
        t = ""
        if a.db:
            ts = [tier(n, spdb, gendb) for n in sorted({name} | d["aliases"])]
            t = min(ts, key=["A", "A1", "B", "C"].index)
        # action: what would add MAT frequency information for this species
        if decision.startswith("NOT_TYPABLE"):
            action = "not_typable"           # homothallic or silent cassettes (causes.tsv): more runs add nothing until active-MAT typing exists (#322)
        elif decision.startswith("PROVISIONAL_HOMOTHALLIC"):
            action = "homothallic"
        elif decision.startswith(PASSING):
            action = "expand" if untyped or loc_new else "complete"
        elif decision.startswith("UNVALIDATED"):
            action = "unvalidated"
        elif decision.startswith("HOLD"):
            action = "blocked_panel" if PANEL_HOLD.search(decision) else "held"
        elif tp:
            action = "typed_no_gate"         # typed with a panel outside the waves (pilot4: Pyricularia oryzae) but never gated; not a reference gap
        elif len(elig) + loc_new < a.min_isolates:
            action = "too_few"
        elif t in ("A", "A1", "B"):
            action = "new_candidate"
        else:
            action = "needs_reference" if a.db else "new_candidate"
        if action == "too_few" and not (tp or tl or loc):
            continue
        largest = max(per_study.values()) if per_study else 0
        ss = d["runs"]
        rows.append({
            "species": name, "species_taxid": "" if key.startswith("name:") else key, "reference_tier": t,
            "wave": gate[2], "gate": decision, "action": action,
            "public_isolates_eligible": len(elig), "public_studies": len(per_study),
            "largest_study_frac": f"{largest / len(elig):.2f}" if elig else "",
            "excluded_isolates": ";".join(f"{k}:{v}" for k, v in sorted(excl.items())),
            "typed_public": len(set(tp)), "typed_local": len(tl), "typed_local_matched_public": len(tl_bs),
            "untyped_public_eligible": len(untyped),
            "local_with_files": len(loc_files), "local_likely_public": len(loc_public), "local_not_public": loc_new,
            "calls": ";".join(f"{c or 'blank'}:{n}" for c, n in sorted(calls.items())),
            "minor_idiomorph_frac": f"{minor / n_single:.3f}" if n_single and len(single) >= 2 else "",
            "ci95_halfwidth_now": wilson_halfwidth(minor, n_single) if len(single) >= 2 else "",
            "ci95_halfwidth_all": wilson_halfwidth(round(minor / n_single * n_all), n_all) if n_single and len(single) >= 2 else "",
            "genome_mb": gsize.get(key, ""),
        })
        if action in ("expand", "unvalidated", "held", "typed_no_gate"):
            slug = gate[3] or re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")     # the wave's own group slug, so a run list finds the panel of plan/<slug>/
            for bs in sorted(untyped):
                r = elig[bs]
                cov = f"{r['bases'] / (gsize[key] * 1e6):.1f}" if key in gsize else ""
                untyped_rows.append([slug, name, r["run_accession"], r["read_count"], r["base_count"], r["library_layout"],
                                     (r.get("strain") or r.get("isolate") or ""), r["study_accession"], bs, cov,
                                     "yes" if MIN_SEL <= r["bases"] <= MAX_SEL else "no", action, gate[2]])
    order = ["expand", "unvalidated", "held", "typed_no_gate", "blocked_panel", "new_candidate", "needs_reference", "homothallic", "not_typable", "complete", "too_few"]
    rows.sort(key=lambda r: (order.index(r["action"]), -r["untyped_public_eligible"] - r["local_not_public"], r["species"]))
    cols = list(rows[0].keys())
    with open(os.path.join(out, "species_availability.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    untyped_rows.sort()
    with gzip.GzipFile(os.path.join(out, "untyped_runs.tsv.gz"), "wb", mtime=0) as gz:
        buf = ["\t".join(["group", "species", "run_accession", "read_count", "base_count", "library_layout", "strain", "study_accession",
                          "biosample", "est_coverage", "in_select_window", "action", "wave"])]
        buf += ["\t".join(map(str, r)) for r in untyped_rows]
        gz.write(("\n".join(buf) + "\n").encode())
    if local_rows:
        with open(os.path.join(out, "local_samples.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(local_rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(sorted(local_rows, key=lambda r: (r["subproject"], r["sample"])))

    inputs = {"runs": a.runs, "taxmap": a.taxmap, "genome_sizes": a.genome_sizes or ""}
    inputs.update({f"typed[{i}]": p for i, p in enumerate(a.typed)})
    inputs.update({f"gates[{i}]": p for i, p in enumerate(a.gates)})
    inputs.update({f"typed_local[{i}]": p for i, p in enumerate(a.typed_local)})
    snap = {"run_table_date": date, "popseq_commit": a.popseq_commit, "min_isolates": a.min_isolates, "skip_species": sorted(skip),
            "min_bases": MIN_BASES, "counts": dict(sorted(counts.items())), "typed_runs_not_in_run_table": dict(sorted(typed_unplaced.items())),
            "local_inventory": bool(a.local_inventory),
            "inputs": {k: {"path": os.path.basename(v), "sha256": sha256(v)} for k, v in sorted(inputs.items()) if v}}
    with open(os.path.join(out, "snapshot.json"), "w") as fh:
        json.dump(snap, fh, indent=1, sort_keys=True)
        fh.write("\n")
    write_summary(out, rows, snap, untyped_rows)
    if a.previous:
        write_changes(out, rows, a.previous)
    print(f"{out}: {len(rows)} species, {len(untyped_rows)} untyped eligible runs; actions:",
          dict(collections.Counter(r["action"] for r in rows)))


def write_summary(out, rows, snap, untyped_rows):
    by = collections.defaultdict(list)
    for r in rows:
        by[r["action"]].append(r)
    meaning = {
        "expand": "passed the wave gate; more public or local isolates exist than were typed (the 50-run cap)",
        "unvalidated": "typed but the gate could not validate the panel; untyped isolates listed, run only after a curator ruling",
        "held": "held for typing-rate or agreement reasons; untyped isolates listed for a re-test",
        "typed_no_gate": "typed with a panel outside the waves (pilot runs) but never gated; untyped isolates listed",
        "blocked_panel": "held because no usable panel could be built; needs reference or panel work, not more reads",
        "new_candidate": "not yet in a wave; MATPredict has records for the species or its genus",
        "needs_reference": "not yet in a wave; no MATPredict record in the genus",
        "homothallic": "passed as homothallic: an idiomorph frequency is not meaningful",
        "not_typable": "not typable by presence (homothallic or silent cassettes, causes.tsv): waits on active-MAT typing (#322)",
        "complete": "passed; every eligible public isolate is typed",
        "too_few": "fewer eligible isolates than the minimum, kept because MATPredict has typed or holds some",
    }
    L = [f"# Reads availability snapshot {snap['run_table_date']}", "",
         f"PopSeq kingdom run table of {snap['run_table_date']} (PopSeq commit `{snap['popseq_commit'] or 'unknown'}`) joined to MATPredict's typed runs. "
         f"Generated by `reads_availability.py`; inputs and their sha256 in `snapshot.json`. Isolate = BioSample; eligible = natural-isolate study "
         f"(survey rule) and >= {snap['min_bases'] // 10**6} Mb.", "",
         "| action | species | eligible public isolates | typed | untyped eligible | local, not public | meaning |", "|---|---|---|---|---|---|---|"]
    for k, v in by.items():
        L.append(f"| `{k}` | {len(v)} | {sum(r['public_isolates_eligible'] for r in v):,} | {sum(r['typed_public'] + r['typed_local'] for r in v):,} | "
                 f"{sum(r['untyped_public_eligible'] for r in v):,} | {sum(r['local_not_public'] for r in v):,} | {meaning[k]} |")
    L += ["", f"`untyped_runs.tsv.gz`: {len(untyped_rows):,} runs, one per untyped eligible isolate, for `expand`, `unvalidated`, `held` and `typed_no_gate` species "
          "(the run nearest 3 Gb inside 0.6-12 Gb, as `select_runs.py`). Columns fit the run list of `CONTRACT.md` once a panel is named.", ""]
    for k, title, n in (("expand", "Largest gains for passed species", 25), ("new_candidate", "Largest new candidates (records for the species or genus)", 25),
                        ("needs_reference", "Largest species with no MATPredict record in the genus", 15), ("unvalidated", "Unvalidated species with data waiting", 15)):
        v = by.get(k, [])
        if not v:
            continue
        L += [f"## {title}", "", "| species | gate | typed | untyped eligible | local, not public | studies | largest study | calls | CI95 half-width now -> all |",
              "|---|---|---|---|---|---|---|---|---|"]
        for r in v[:n]:
            ci = f"{r['ci95_halfwidth_now']} -> {r['ci95_halfwidth_all']}" if r["ci95_halfwidth_now"] else ""
            L.append(f"| *{r['species']}* | {r['gate']} | {r['typed_public'] + r['typed_local']} | {r['untyped_public_eligible']:,} | {r['local_not_public']} | "
                     f"{r['public_studies']} | {r['largest_study_frac']} | {r['calls']} | {ci} |")
        L.append("")
    L += ["## Read with care", "",
          "- `largest_study_frac` near 1 means one study supplies most isolates; outbreak and surveillance sets are often clonal, so their idiomorph counts are not a population frequency.",
          "- The CI columns are binomial (Wilson) half-widths for the minor idiomorph, now and if every untyped eligible isolate were typed at the species' current typing rate; they ignore clonality.",
          "- `local_likely_public` matches a local sample to the archive by accession, BioSample or strain name. A strain-name match is weak (spec 005 section 4); `local_not_public` is an upper bound on new strains.",
          "- Coccidioides strains typed from local files are matched to public isolates by strain name only.",
          "- Species are keyed by NCBI species taxid (PopSeq's map), so synonyms in the archive (var., f. sp., renamed genera) are pooled; the survey of 2026-10-26 pooled by the first two words of the name.", ""]
    open(os.path.join(out, "SUMMARY.md"), "w").write("\n".join(L))


def write_changes(out, rows, prev_dir):
    p = os.path.join(prev_dir, "species_availability.tsv")
    if not os.path.exists(p):
        print("no species_availability.tsv in", prev_dir, "- CHANGES.md not written")
        return
    prev = {(r["species_taxid"] or r["species"]): r for r in read_tsv(p)}
    cur = {(r["species_taxid"] or r["species"]): r for r in rows}
    new = [cur[k] for k in cur if k not in prev and cur[k]["action"] != "too_few"]
    gone = [prev[k] for k in prev if k not in cur]
    moved = [(prev[k], cur[k]) for k in cur if k in prev and prev[k]["action"] != cur[k]["action"]]
    grew = sorted(((cur[k]["public_isolates_eligible"] - int(prev[k]["public_isolates_eligible"]), cur[k]) for k in cur if k in prev),
                  key=lambda x: (-x[0], x[1]["species"]))
    grew = [g for g in grew if g[0] > 0]
    shrank = [(int(prev[k]["public_isolates_eligible"]) - cur[k]["public_isolates_eligible"], cur[k]) for k in cur if k in prev
              and cur[k]["public_isolates_eligible"] < int(prev[k]["public_isolates_eligible"])]
    L = [f"# Changes since {os.path.basename(os.path.normpath(prev_dir))}", "",
         f"- {len(new)} species new to the table, {len(gone)} gone, {len(moved)} changed action, {len(grew)} gained eligible isolates, {len(shrank)} lost some.", ""]
    if shrank:
        L += ["## Lost eligible isolates (check: archives rarely drop records; usually the query or a filter changed)", ""]
        L += [f"- *{r['species']}*: -{n}" for n, r in sorted(shrank, key=lambda x: -x[0])[:40]] + [""]
    if moved:
        L += ["## Changed action", ""] + [f"- *{c['species']}*: {p['action']} -> {c['action']}" for p, c in moved] + [""]
    if grew:
        L += ["## Gained eligible isolates", ""] + [f"- *{r['species']}* ({r['action']}): +{n}" for n, r in grew[:60]] + [""]
    if new:
        L += ["## New species", ""] + [f"- *{r['species']}* ({r['action']}): {r['public_isolates_eligible']} eligible" for r in
                                       sorted(new, key=lambda r: -r["public_isolates_eligible"])[:60]] + [""]
    open(os.path.join(out, "CHANGES.md"), "w").write("\n".join(L))


if __name__ == "__main__":
    main()
