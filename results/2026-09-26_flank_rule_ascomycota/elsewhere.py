"""For calls whose core best hits are elsewhere: does that best hit fall inside ANOTHER
detected call (core-modelled or not) or a suppressed locus in the same report?"""
import csv, collections, subprocess, yaml
from yaml import CSafeLoader as L
R = list(csv.DictReader(open("classified.tsv"), delimiter="\t"))
base = {"cap6": "../2026-09-26_polish_cap/cap6/runs"}
MOD = {"polished_agree", "polished_disagree", "polished_single"}
out = csv.writer(open("elsewhere.tsv", "w"), delimiter="\t")
out.writerow(["genome", "order", "outcome", "verdict", "call", "best_hit_locations"])
cnt = collections.Counter()
for r in R:
    rd = base.get(r["panel"], f"../2026-09-26_early_diverging/{r['panel']}/runs") + "/" + r["genome"]
    rep = yaml.load(open(rd + "/detection_report.yaml"), Loader=L)
    loci = []
    for d in rep.get("detected") or []:
        if d["contig"] == r["contig"] and str(d["start"]) == r["start"] and d["family"] == r["mat_family"]: continue
        cm = any(e["role"] == "core_MAT" and (e.get("status") in MOD or e.get("method") == "diamond_proteome") for e in d.get("gene_evidence") or [])
        loci.append(("called_coremodelled" if cm else "called_flankcarried", d["contig"], d["start"], d["end"]))
    for d in rep.get("suppressed_loci") or []:
        loci.append(("withheld_by_bar", d["contig"], d["start"], d["end"]))
    txt = subprocess.run(["zstdcat", f"hits/{r['genome']}.tsv.zst"], capture_output=True, text=True).stdout
    best = {}
    for l in txt.splitlines():
        q, s, pid, ln, qs, qe, ss, se, ev, bs, ql = l.split("\t")
        gname = q.split("|")[-1]; bs = float(bs)
        if float(ev) > 1e-5: continue
        if gname not in best or bs > best[gname][0]: best[gname] = (bs, s, *sorted((int(ss), int(se))), float(pid))
    locs = []
    for gname in r["genes_best_elsewhere"].split(",") if r["genes_best_elsewhere"] else []:
        bs, s, lo, hi, pid = best[gname]
        where = [k for k, c, a, b in loci if c == s and a - 20000 <= lo and hi <= b + 20000]
        locs.append(f"{gname}@{s}:{lo}({pid:.0f}%)={'|'.join(sorted(set(where))) or 'no_locus'}")
    if r["verdict"] == "best_hit_elsewhere":
        tag = "coremodelled_call" if any("called_coremodelled" in x for x in locs) else ("other_locus" if any("=no_locus" not in x for x in locs) else "no_locus")
        cnt[(r["order"], r["outcome"], tag)] += 1
    out.writerow([r["genome"], r["order"], r["outcome"], r["verdict"], f"{r['contig']}:{r['start']}-{r['end']}", ";".join(locs)])
for k, v in sorted(cnt.items()): print(v, k)
