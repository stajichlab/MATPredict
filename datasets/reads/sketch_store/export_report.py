#!/usr/bin/env python3
"""Population MAT report of a sketch store, for the site (web/intro/reads/population.html): one JSON and one TSV of public runs.

    export_report.py --store STORE --matfreq MATFREQ_DIR --species "Aspergillus fumigatus" --slug aspergillus_fumigatus \
        [--titles datasets/reads/sra_survey/pilot/study_titles.tsv] --out datasets/reads/sketch_store/reports/aspergillus_fumigatus

Runs on the HPCC where the store lives; its output is small and committed, and scripts/build_reads_site.py copies it into the site.
Inputs: STORE/manifest.tsv, STORE/query/{mat,qc,distances,screen_*}, STORE/exclude.tsv, and the output folder of mat_frequency.py.
Local (Population_Genomics) samples count in every number but are never listed by name: the run table holds public runs only
and the PCoA points carry no identifiers. Writes OUT/population.json and OUT/runs_public.tsv."""
import argparse, collections, csv, datetime, glob, json, os, re, statistics

csv.field_size_limit(10 ** 9)


def rd(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if path and os.path.exists(path) else []


def num(v, d=4):
    try:
        return round(float(v), d)
    except (TypeError, ValueError):
        return None


def freq(path):
    out = []
    for r in rd(path):
        out.append(dict(group=r["group"], n=int(r["n_typed"]), MAT1_1=int(r["MAT1_1"]), MAT1_2=int(r["MAT1_2"]), both=int(r["both"]),
                        frac=num(r["frac_MAT1_1"], 3), lo=num(r["ci95_lo"], 3), hi=num(r["ci95_hi"], 3), p=num(r["p_vs_0_5"], 4)))
    return out


def screens(qdir):
    """Target screens: per target the store-wide present/partial/absent counts and the comparison with a truth set, from summary.txt."""
    rows = []
    for f in sorted(glob.glob(os.path.join(qdir, "screen_*", "summary.txt"))):
        for line in open(f):
            m = re.match(r"(\S+)\s+len\s+(\d+)\s+specific hashes\s+(\d+)\s+present (\d+)\s+partial (\d+)\s+absent (\d+)(?:\s+\|\s+vs truth TP (\d+) FN (\d+) FP (\d+) TN (\d+))?", line)
            if m:
                g = [int(x) if x else None for x in m.groups()[1:]]
                rows.append(dict(screen=os.path.basename(os.path.dirname(f))[7:], target=m.group(1), length=g[0], specific_hashes=g[1], present=g[2],
                                 partial=g[3], absent=g[4], truth=None if g[5] is None else dict(TP=g[5], FN=g[6], FP=g[7], TN=g[8])))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--store", required=True); ap.add_argument("--matfreq", required=True); ap.add_argument("--species", required=True)
    ap.add_argument("--slug", required=True); ap.add_argument("--titles"); ap.add_argument("--out", required=True)
    ap.add_argument("--species-id", help="species_id.tsv (sample_id, best, best_ani): nearest reference genome of the isolates unlike the rest")
    ap.add_argument("--clone", default="0.9995"); ap.add_argument("--lineage", default="0.998"); ap.add_argument("--min-group", type=int, default=10)
    a = ap.parse_args(argv)
    q = os.path.join(a.store, "query")
    man = {r["sample_id"]: r for r in rd(os.path.join(a.store, "manifest.tsv"))}
    mat = {r["sample_id"]: r for r in rd(f"{q}/mat/mat_calls.tsv")}
    sq = {r["sample_id"]: r for r in rd(f"{q}/distances/samples_qc.tsv")}
    qc = {r["sample_id"]: r for r in rd(f"{q}/qc/qc.tsv")}
    fs = {r["sample_id"]: r for r in rd(os.path.join(a.matfreq, "samples.tsv"))}
    titles = {r["study"]: r["title"] for r in rd(a.titles)} if a.titles else {}
    os.makedirs(a.out, exist_ok=True)

    src = collections.Counter(r["source"] for r in man.values())
    sketched = sum(1 for s in man if os.path.exists(os.path.join(a.store, "sigs", s + ".sig.zip")))
    excluded = collections.Counter(r["excluded"] or "kept" for r in fs.values())
    overall = {k: freq(os.path.join(a.matfreq, f"freq_overall_{k}.tsv"))[0] for k in ("samples", "clone-corrected")}
    conc = collections.Counter((m.get("reads_type_call", ""), m.get("sketch_call", "")) for m in mat.values())

    # lineages that depart from 0.5 (clone-corrected, p < 0.01): what they are made of
    lin_rows = freq(os.path.join(a.matfreq, "freq_lineage_clones.tsv"))
    kept = [r for r in fs.values() if not r["excluded"]]
    for L in lin_rows:
        mem = [r for r in kept if (r["lineage"] or "none") == L["group"]]
        L["samples"] = len(mem); L["clones"] = len({r["clone"] for r in mem})
        L["studies"] = collections.Counter(r["study"] for r in mem).most_common(4)
        L["countries"] = collections.Counter(r["country"] for r in mem).most_common(4)

    # distance outliers in four kinds: too little data (mean ANI < 0.95 with few k-mers: a few hundred reads), another species (mean ANI
    # < 0.95 otherwise), poor sample (a QC flag or no reads-type depth), divergent lineage
    ani = [float(r["mean_ani"]) for r in sq.values() if r.get("mean_ani")]
    outl = [s for s, r in sq.items() if r.get("outlier") == "1"]
    kind = {}
    for s in outl:
        ex = fs.get(s, {}).get("excluded", "")
        if ex in ("too_little_data", "other_species"):
            kind[s] = ex
        elif float(sq[s]["mean_ani"]) < 0.95:
            kind[s] = "other_species"      # a curator exclusion of another species (local Moose_3-1.1 is A. nidulans)
        elif (qc.get(s, {}).get("flag") or "") or mat.get(s, {}).get("reads_type_call") in ("low_depth", "none"):
            kind[s] = "poor_sample"
        else:
            kind[s] = "divergent_lineage"
    spid = {r["sample_id"]: r for r in rd(a.species_id)} if a.species_id else {}

    def by_study(reason):
        g = collections.defaultdict(list)
        for s, f in fs.items():
            if f["excluded"] == reason or (reason == "other_species" and f["excluded"] == "curator_exclude" and spid.get(s, {}).get("best", "").find("fumigatus") < 0 and s in spid):
                g[man.get(s, {}).get("study_accession", "")].append(s)
        rows = []
        for st, ss in sorted(g.items(), key=lambda x: -len(x[1])):
            best = collections.defaultdict(list)      # nearest reference -> implied identities (one entry per isolate)
            for s in ss:
                if s in spid:
                    v = float(spid[s].get("best_ani") or 0)
                    best[spid[s]["best"] if v >= 0.95 else "unresolved"].append(v)
            label = lambda k, v: k if k == "unresolved" else f"{k}, ANI {min(v):.3f}" + (f"-{max(v):.3f}" if max(v) > min(v) else "")
            best = collections.Counter({label(k, v) + (f" x{len(v)}" if len(v) > 1 else ""): len(v) for k, v in best.items()})
            rows.append(dict(study=st if man_public(st, man) else "local", n=len(ss), mean_ani=round(statistics.mean(float(fs[s]["mean_ani"]) for s in ss if fs[s]["mean_ani"]), 3),
                             title=titles.get(st, ""), nearest="; ".join(k for k, _ in best.most_common()),
                             reads=sum(int(man[s].get("read_count") or 0) for s in ss) // max(1, len(ss))))
        return rows
    other_studies, thin_studies = by_study("other_species"), by_study("too_little_data")

    # PCoA points: PC1, PC2, call, public/local; no identifiers
    var = None
    for lg in sorted(glob.glob(os.path.join(a.store, "logs", "query.*.log")), key=os.path.getmtime, reverse=True):
        m = re.search(r"variance share of axes 1-3 ([\d.]+), ([\d.]+), ([\d.]+)", open(lg).read())
        if m:
            var = [float(x) for x in m.groups()]; break
    code = {"MAT1-1": 1, "MAT1-2": 2, "both": 3}
    pts = []
    for r in rd(f"{q}/distances/pcoa.tsv"):
        s = r["sample_id"]; f = fs.get(s, {})
        if f.get("excluded") in ("other_species", "too_little_data", "curator_exclude"):
            continue
        pts.append([round(float(r["PCo1"]), 6), round(float(r["PCo2"]), 6), code.get(f.get("call"), 0), 1 if man.get(s, {}).get("source") == "local" else 0])

    report = dict(
        slug=a.slug, species=a.species, built=datetime.date.today().isoformat(),
        store=dict(samples=len(man), public=src["public"], local=src["local"], sketched=sketched, missing=len(man) - sketched,
                   params="k=21,k=31 scaled=100 abund, k-mers seen once dropped; 30x subsample", clone_ani=float(a.clone), lineage_ani=float(a.lineage)),
        excluded=dict(excluded),
        overall=overall,
        concordance=[dict(reads_type=k[0], sketch=k[1], n=v) for k, v in conc.most_common()],
        by_country=dict(samples=freq(os.path.join(a.matfreq, "freq_country_samples.tsv")), clones=freq(os.path.join(a.matfreq, "freq_country_clones.tsv"))),
        by_study=dict(samples=freq(os.path.join(a.matfreq, "freq_study_samples.tsv")), clones=freq(os.path.join(a.matfreq, "freq_study_clones.tsv"))),
        by_source=dict(samples=freq(os.path.join(a.matfreq, "freq_source_samples.tsv")), clones=freq(os.path.join(a.matfreq, "freq_source_clones.tsv"))),
        lineages=lin_rows, min_group=a.min_group,
        outliers=dict(n=len(outl), median_ani=round(statistics.median(ani), 5) if ani else None,
                      cutoff=round(max(float(sq[s]["mean_ani"]) for s in outl), 5) if outl else None,
                      kinds=dict(collections.Counter(kind.values())), other_species_studies=other_studies, too_little_data_studies=thin_studies),
        pcoa=dict(var=var, points=pts),
        screens=screens(q),
        study_titles={r["group"]: titles[r["group"]] for r in freq(os.path.join(a.matfreq, "freq_study_samples.tsv")) if r["group"] in titles},
    )
    json.dump(report, open(os.path.join(a.out, "population.json"), "w"), separators=(",", ":"))
    cols = ["biosample", "run", "strain", "study", "country", "first_public", "call", "sketch_call", "reads_type_call", "clone", "lineage", "mean_ani", "excluded"]
    with open(os.path.join(a.out, "runs_public.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(cols)
        for s, f in sorted(fs.items()):
            m = man.get(s, {})
            if m.get("source") != "public":
                continue
            w.writerow([s, m.get("run_accession", ""), m.get("strain", ""), m.get("study_accession", ""), f["country"], m.get("first_public", ""), f["call"],
                        mat.get(s, {}).get("sketch_call", ""), mat.get(s, {}).get("reads_type_call", ""), f["clone"], f["lineage"], f["mean_ani"], f["excluded"]])
    print(f"{a.slug}: {len(fs)} samples, {excluded.get('kept', 0)} kept; clone-corrected MAT1-1 {overall['clone-corrected']['frac']}; "
          f"{len(pts)} PCoA points; {len(report['screens'])} screen targets -> {a.out}")


def man_public(study, man):
    """True when the study is a public accession (local samples' 'study' is a folder name of the Population_Genomics project)."""
    return bool(re.match(r"^(PRJ|[SED]RP)", study or ""))


if __name__ == "__main__":
    main()
