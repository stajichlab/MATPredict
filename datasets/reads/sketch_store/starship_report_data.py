#!/usr/bin/env python3
"""Collect the outputs of starship_dynamics.py (and the screen runs they came from) into one JSON for the Starship report page.

    starship_report_data.py --default DIR [--tuned DIR] --candidates TSV --before SUMMARY --after SUMMARY [--dropped TSV]
                            [--baf-family JSON] [--meta key=value ...] --out report_data.json [--template report.html --html OUT.html]

Every table is read from a named file, nothing is recomputed except the 2x2 bafA / Nebuchadnezzar-h1 table and the heatmap order, so a
number on the page can be traced to a file in DIR. --meta adds provenance strings (commit, job ids, paths) shown in the page footer.
With --template/--html the JSON is written into the template at the marker /*__DATA__*/null."""
import argparse, csv, json, os, re, sys

LINE = re.compile(r"^(\S+)\s+len\s+(\d+)\s+specific hashes\s+(\d+)\s+present\s+(\d+)\s+partial\s+(\d+)\s+absent\s+(\d+)"
                  r"(?:\s+\|\s+vs truth TP (\d+) FN (\d+) FP (\d+) TN (\d+))?(?:\s+\|\s+dropped (.*))?$")


def tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        try:
            return float(x)
        except (TypeError, ValueError):
            return x


def conv(rows):
    return [{k: num(v) for k, v in r.items()} for r in rows]


def screen_summary(path):
    out = {}
    for line in open(path):
        m = LINE.match(line.strip())
        if m:
            g = m.groups()
            out[g[0]] = dict(len=int(g[1]), specific=int(g[2]), present=int(g[3]), partial=int(g[4]), absent=int(g[5]),
                             truth=dict(zip(("TP", "FN", "FP", "TN"), map(int, g[6:10]))) if g[6] else None, dropped=g[10] or "")
    return out


def collect(d):
    out = dict(prevalence=conv(tsv(f"{d}/prevalence.tsv")), parsimony=conv(tsv(f"{d}/parsimony.tsv")),
               by_lineage=conv(tsv(f"{d}/prevalence_by_lineage.tsv")), by_country=conv(tsv(f"{d}/prevalence_by_country.tsv")))
    out["summary"] = open(f"{d}/summary.txt").read().splitlines()
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--default", required=True); ap.add_argument("--tuned"); ap.add_argument("--candidates", required=True)
    ap.add_argument("--before", required=True); ap.add_argument("--after", required=True); ap.add_argument("--dropped")
    ap.add_argument("--meta", nargs="*", default=[]); ap.add_argument("--out", required=True)
    ap.add_argument("--baf-family", help="baf_family_summary.json from baf_family.py summarize (adds the baf family section)")
    ap.add_argument("--pair", default="bafA,Nebuchadnezzar-h1", help="two elements for the clone-level 2x2 table")
    ap.add_argument("--template"); ap.add_argument("--html")
    a = ap.parse_args(argv)
    cm = tsv(f"{a.default}/clone_matrix.tsv")
    fixed = ("clone", "n_members", "representative", "strain", "study", "country", "lineage")
    elements = [c for c in cm[0] if c not in fixed]
    rank = {r["clone"]: int(r["rank"]) for r in tsv(f"{a.default}/tree_order.tsv")} if os.path.exists(f"{a.default}/tree_order.tsv") else {}
    cm.sort(key=lambda r: (rank.get(r["clone"], 10 ** 9), r["clone"]))
    code = {"1": "1", "0": "0", "NA": "2"}
    heat = [dict(clone=r["clone"], n=int(r["n_members"]), lineage=r["lineage"], country=r["country"], in_tree=r["clone"] in rank,
                 bits="".join(code[r[e]] for e in elements)) for r in cm]
    # the two --pair elements, clone level (typed in both): keys are "<first><second>" presence, e.g. "11" = both
    e1, e2 = a.pair.split(",")
    two = None
    if e1 in elements and e2 in elements:
        two = {"11": 0, "10": 0, "01": 0, "00": 0}
        for r in cm:
            if r[e1] != "NA" and r[e2] != "NA":
                two[r[e1] + r[e2]] += 1
    sweep = {}
    if os.path.exists(f"{a.default}/threshold_sweep.tsv"):
        for r in conv(tsv(f"{a.default}/threshold_sweep.tsv")):
            sweep.setdefault(r["element"], []).append({k: r[k] for k in ("present_cut", "TP", "FN", "FP", "TN", "precision", "recall")})
    data = dict(
        meta=dict(kv.split("=", 1) for kv in a.meta), elements=elements, pair=[e1, e2], default=collect(a.default),
        tuned=collect(a.tuned) if a.tuned else None, heat=heat, baf_neb=two, sweep=sweep,
        cooccurrence=conv(tsv(f"{a.default}/cooccurrence.tsv")), candidates=conv(tsv(a.candidates)),
        baf_family=json.load(open(a.baf_family)) if a.baf_family else None,
        rescreen=dict(before=screen_summary(a.before), after=screen_summary(a.after), dropped=conv(tsv(a.dropped)) if a.dropped else []))
    with open(a.out, "w") as fh:
        json.dump(data, fh, separators=(",", ":"))
    if a.template and a.html:
        t = open(a.template).read()
        if "/*__DATA__*/null" not in t:
            sys.exit("template has no /*__DATA__*/null marker")
        open(a.html, "w").write(t.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")).replace("</", "<\\/"), 1))
    print(f"{len(elements)} elements, {len(heat)} clone rows ({sum(h['in_tree'] for h in heat)} in the tree), "
          f"{len(data['candidates'])} candidates, bafA x Neb-h1 {two}")


if __name__ == "__main__":
    sys.exit(main())
