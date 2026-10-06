"""Held-out rerun on PR #9 52b3ff9: per-genome tables, by-species summaries,
LCG file-label scoring and the Zygo 23 subset. Read-only analysis.

Old runs (076afe4): ../2026-09-28_lcg_holdout/runs2, ../2026-09-28_mucor_jena_holdout/runs_{scaffolds,contigs}
New runs (52b3ff9): lcg_runs, jena_scaffolds, jena_contigs
"""
import collections, csv, glob, os, re
import yaml

R = os.path.dirname(os.path.abspath(__file__))
P = os.path.dirname(R)
OLD = {"lcg": f"{P}/2026-09-28_lcg_holdout/runs2",
       "jena_scaffolds": f"{P}/2026-09-28_mucor_jena_holdout/runs_scaffolds",
       "jena_contigs": f"{P}/2026-09-28_mucor_jena_holdout/runs_contigs"}
NEW = {"lcg": f"{R}/lcg_runs", "jena_scaffolds": f"{R}/jena_scaffolds",
       "jena_contigs": f"{R}/jena_contigs"}
CONF = {"high": 3, "medium": 2, "low": 1}
# Curator ruling 2026-10-01 (B8): misidentified genomes are listed in the
# versioned db/taxon_overrides.tsv. A genome is excluded from species scoring
# if it is there OR its curator note says MISIDENTIFIED.
OVERRIDES = os.environ.get(
    "MATPREDICT_TAXON_OVERRIDES",
    "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db/taxon_overrides.tsv")


def load_overrides(path):
    if not os.path.exists(path):
        raise SystemExit(f"taxon overrides not found: {path}")
    rows = [l for l in open(path) if l.strip() and not l.startswith("#")]
    return {r["genome_id"]: r for r in csv.DictReader(rows, delimiter="\t")}


def load(d):
    out = {}
    for f in glob.glob(f"{d}/*/detection_report.yaml"):
        g = f.split("/")[-2]
        try:
            out[g] = yaml.safe_load(open(f)) or {}
        except Exception as e:  # noqa: BLE001
            out[g] = {"_error": str(e)}
    return out


def summ(r):
    """Genome-level summary of one report."""
    if r is None:
        return dict(status="no_report")
    dets = [d for d in (r.get("detected") or []) if d.get("family", "").startswith("Mucoromycota")]
    labs = sorted({d.get("idiomorph") for d in dets})
    det = [l for l in labs if l in ("Plus", "Minus")]
    if not dets:
        call = "uncalled"
    elif set(det) == {"Plus", "Minus"}:
        call = "both"
    elif det:
        call = det[0]
    else:
        call = "undetermined"
    best = max((CONF.get(d.get("confidence"), 0) for d in dets), default=0)
    bestconf = {3: "high", 2: "medium", 1: "low", 0: ""}[best]
    wr = collections.Counter(s.get("withheld_reason") for s in (r.get("suppressed_loci") or []))
    ti = r.get("two_idiomorphs") or []
    arr = ti[0].get("arrangement") if ti else ""
    causes = ";".join(ti[0].get("supported_causes") or []) if ti else ""
    calls = ";".join(f"{d.get('contig')}:{d.get('start')}-{d.get('end')}|{d.get('idiomorph')}|"
                     f"{d.get('confidence')}|{d.get('locus_class')}|"
                     f"{(d.get('idiomorph_classifier') or {}).get('classifier_input','')}|"
                     f"{(d.get('idiomorph_classifier') or {}).get('margin','')}"
                     for d in dets)
    return dict(status="called" if dets else "uncalled", call=call, n_calls=len(dets),
                best_conf=bestconf, calls=calls, two_idiomorphs=arr, two_idio_causes=causes,
                gate=wr.get("mat_gene_gate", 0), paralog=wr.get("paralog_class", 0),
                floor=wr.get("below_fraction_floor", 0),
                verification=";".join(sorted({d.get("verification") or "" for d in dets} - {""})))


def file_species(org):
    t = org.split("_")
    sp = [t[0]]
    for x in t[1:]:
        if re.match(r"^(NRRL|RSA|CBS|IMI|URM|ATCC|BCRC|UCR|\d)", x):
            break
        sp.append(x)
        if len(sp) >= 2 and sp[-1] not in ("var.", "f.", "sp.", "aff.", "cf."):
            if len(sp) >= 2 and sp[-2] not in ("var.", "f.", "aff.", "cf."):
                break
    return " ".join(sp)


def write(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def species_summary(rows, key, out):
    by = collections.defaultdict(list)
    for r in rows:
        if r.get("exclude_species"):
            continue
        by[r[key]].append(r)
    res = []
    for sp, rs in sorted(by.items()):
        c = collections.Counter(r["new_call"] for r in rs)
        conf = collections.Counter(r["new_best_conf"] for r in rs if r["new_best_conf"])
        res.append(dict(species=sp, genomes=len(rs), called=sum(r["new_status"] == "called" for r in rs),
                        Plus=c["Plus"], Minus=c["Minus"], both=c["both"],
                        undetermined=c["undetermined"], uncalled=c["uncalled"],
                        high=conf["high"], medium=conf["medium"], low=conf["low"],
                        two_idiomorphs=";".join(sorted({r["new_two_idiomorphs"] for r in rs} - {""})),
                        gate_withheld=sum(int(r["new_gate"]) for r in rs),
                        paralog_withheld=sum(int(r["new_paralog"]) for r in rs),
                        floor_withheld=sum(int(r["new_floor"]) for r in rs)))
    write(out, res, list(res[0].keys()))
    return res


def main():
    old = {k: load(v) for k, v in OLD.items()}
    new = {k: load(v) for k, v in NEW.items()}
    lines = []
    # ---------------- LCG
    global TAXON_OV
    TAXON_OV = load_overrides(OVERRIDES)
    ct = {r["org"]: r for r in csv.DictReader(open(f"{P}/2026-09-28_lcg_holdout/curator_table.tsv"), delimiter="\t")}
    pg = {r["org"]: r for r in csv.DictReader(open(f"{P}/2026-09-28_lcg_holdout/per_genome.tsv"), delimiter="\t")}
    pil = {}
    pf = f"{P}/2026-09-28_lcg_holdout/protein_identity_leak.tsv"
    if os.path.exists(pf):
        for r in csv.DictReader(open(pf), delimiter="\t"):
            pil[r.get("org") or list(r.values())[0]] = r
    truth = {}
    for l in open(f"{P}/2026-09-23_zygo_bar/zygo_truth.tsv"):
        f = l.rstrip("\n").split("\t")
        truth[f[0]] = (f[1], int(f[2]), int(f[3]), f[4])
    rows = []
    for org in sorted(ct):
        c = ct[org]
        o, n = summ(old["lcg"].get(org)), summ(new["lcg"].get(org))
        mis = "MISIDENTIFIED" in (c.get("curator_notes") or "") or org in TAXON_OV
        sp = c.get("curator_taxonomy") or file_species(org)
        row = dict(org=org, species=sp, species_source="curator" if c.get("curator_taxonomy") else "file_name",
                   leakage=pg.get(org, {}).get("leakage", ""),
                   protein_identity_leak="yes" if org in pil else "",
                   label=c.get("curator_mating_type", ""),
                   label_disputed="yes" if "DISPUTED" in (c.get("curator_mating_source") or "") else "",
                   misidentified="yes" if mis else "", exclude_species="yes" if mis else "",
                   **{f"old_{k}": v for k, v in o.items()}, **{f"new_{k}": v for k, v in n.items()})
        rows.append(row)
    cols = list(rows[0].keys())
    write(f"{R}/lcg_per_genome.tsv", rows, cols)
    sp = species_summary(rows, "species", f"{R}/lcg_by_species.tsv")
    called_old = sum(r["old_status"] == "called" for r in rows)
    called_new = sum(r["new_status"] == "called" for r in rows)
    lines.append(f"LCG genomes {len(rows)}; reports new {sum(r['new_status']!='no_report' for r in rows)}; "
                 f"called old {called_old} -> new {called_new}")
    lines.append("LCG new call split: " + str(dict(collections.Counter(r["new_call"] for r in rows))))
    lines.append("LCG old call split: " + str(dict(collections.Counter(r["old_call"] for r in rows))))
    lines.append("LCG new best confidence: " + str(dict(collections.Counter(r["new_best_conf"] for r in rows if r["new_status"] == "called"))))
    lines.append("LCG two_idiomorphs (new): " + str(dict(collections.Counter(r["new_two_idiomorphs"] for r in rows if r["new_two_idiomorphs"]))))
    lines.append(f"LCG species summarised {len(sp)} (misidentified excluded: {sum(1 for r in rows if r['misidentified'])})")
    # label scoring
    def score(rows_, tag):
        res = collections.Counter()
        det = []
        for r in rows_:
            lab, call = r["label"], r["new_call"]
            if call == "uncalled":
                v = "uncalled"
            elif call == lab:
                v = "agree"
            elif call == "both":
                v = "both"
            elif call == "undetermined":
                v = "undetermined"
            else:
                v = "disagree"
            res[(v, r["new_best_conf"] or "-")] += 1
            res[v] += 1
            det.append((r["org"], lab, call, r["new_best_conf"], v))
        lines.append(f"  {tag}: n={len(rows_)} agree={res['agree']} disagree={res['disagree']} uncalled={res['uncalled']} "
                     f"both={res['both']} undetermined={res['undetermined']}; by tier " +
                     str({k: v for k, v in res.items() if isinstance(k, tuple)}))
        for d in det:
            if d[4] != "agree":
                lines.append(f"    {d[4]}: {d[0]} label={d[1]} call={d[2]} conf={d[3]}")
    lab_rows = [r for r in rows if r["label"]]
    clean = [r for r in lab_rows if not r["label_disputed"] and not r["misidentified"]
             and r["leakage"] not in ("training_leak", "zygo23")]
    lines.append("LCG file-label scoring (new code):")
    score(clean, "clean (no disputed/misidentified/training/zygo)")
    score([r for r in lab_rows if r["leakage"] == "zygo23"], "labelled zygo23")
    score([r for r in lab_rows if r["label_disputed"]], "disputed (for reference)")
    # zygo
    zok_i = zok_l = 0
    zbad = []
    for org, (sc, s, e, idio) in truth.items():
        rep = new["lcg"].get(org) or {}
        dets = [d for d in (rep.get("detected") or []) if d.get("family", "").startswith("Mucoromycota")]
        loc = [d for d in dets if d.get("contig") == sc]
        li = bool(loc)
        ii = any(d.get("idiomorph") == idio for d in loc)
        zok_l += li
        zok_i += ii
        if not (li and ii):
            zbad.append((org, sc, idio, [(d.get("contig"), d.get("idiomorph")) for d in dets]))
    lines.append(f"Zygo 23 subset (LCG scaffolds, new code): locus {zok_l}/{len(truth)}, idiomorph {zok_i}/{len(truth)}")
    for z in zbad:
        lines.append(f"  zygo miss: {z}")
    # ---------------- Jena
    jt = {r["strain"]: r for r in csv.DictReader(open(f"{P}/2026-09-28_mucor_jena_holdout/curator_table.tsv"), delimiter="\t")}
    jrows = []
    for s in sorted(jt):
        c = jt[s]
        row = dict(strain=s, species=c.get("curator_taxonomy") or "(not in curator table)",
                   type_status=c.get("curator_type_status", ""), leakage=c.get("leakage", "") or c.get("training_leak", ""))
        for inp in ("scaffolds", "contigs"):
            o, n = summ(old[f"jena_{inp}"].get(s)), summ(new[f"jena_{inp}"].get(s))
            for k, v in o.items():
                row[f"old_{inp}_{k}"] = v
            for k, v in n.items():
                row[f"new_{inp}_{k}"] = v
        # species summary uses scaffolds
        for k in ("status", "call", "best_conf", "two_idiomorphs", "gate", "paralog", "floor"):
            row[f"new_{k}"] = row.get(f"new_scaffolds_{k}", "")
        jrows.append(row)
    write(f"{R}/jena_per_strain.tsv", jrows, list(jrows[0].keys()))
    jsp = species_summary(jrows, "species", f"{R}/jena_by_species.tsv")
    lines.append(f"Jena strains {len(jrows)}; scaffold reports new {sum(r['new_scaffolds_status']!='no_report' for r in jrows)}; "
                 f"called old {sum(r['old_scaffolds_status']=='called' for r in jrows)} -> new {sum(r['new_scaffolds_status']=='called' for r in jrows)}")
    lines.append("Jena new call split (scaffolds): " + str(dict(collections.Counter(r["new_scaffolds_call"] for r in jrows))))
    lines.append("Jena old call split (scaffolds): " + str(dict(collections.Counter(r["old_scaffolds_call"] for r in jrows))))
    agree_sc = sum(r["new_scaffolds_call"] == r["new_contigs_call"] for r in jrows if r["new_contigs_status"] != "no_report")
    lines.append(f"Jena scaffolds vs contigs same call (new): {agree_sc}/{sum(r['new_contigs_status']!='no_report' for r in jrows)}")
    lines.append("Jena two_idiomorphs (new, scaffolds): " + str([(r['strain'], r['species'], r['new_scaffolds_two_idiomorphs']) for r in jrows if r['new_scaffolds_two_idiomorphs']]))
    lines.append("Jena strains whose scaffold call changed (old -> new):")
    for r in jrows:
        if (r["old_scaffolds_call"], r["old_scaffolds_best_conf"]) != (r["new_scaffolds_call"], r["new_scaffolds_best_conf"]):
            lines.append(f"  {r['strain']} ({r['species']}): {r['old_scaffolds_call']}/{r['old_scaffolds_best_conf']} -> "
                         f"{r['new_scaffolds_call']}/{r['new_scaffolds_best_conf']}  paralog_withheld={r['new_scaffolds_paralog']}")
    for s in ("CBS221_71", "CBS223_63", "CBS763_74"):
        r = next((x for x in jrows if x["strain"] == s), None)
        if r:
            lines.append(f"  P1 case {s} ({r['species']}): new {r['new_scaffolds_calls']} paralog_withheld={r['new_scaffolds_paralog']}")
    lines.append(f"Jena species summarised {len(jsp)}")
    open(f"{R}/summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
