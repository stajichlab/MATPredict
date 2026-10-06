"""Replay candidate rules that would withhold weak sexM-like calls on the gate's
FLANK route (winning classifier score below the build threshold, or fragment-typed).

Rules (applied only to calls that pass the current gate via flank support):
  R1  the supporting flanks must include rnhA
  R2a the winning core model must be >= 100 aa (fragment-typed calls count as 0 aa)
  R2b ... >= 125 aa
  R3  R1 or R2a must hold (withhold only if neither)
Current gate: model score >= THRESH, else >= 2 distinct roster flanks modelled at >= 40%.
Writes replay_rules.tsv (one row per call) and prints per-set summaries.
"""
import collections, csv, glob, os, sys
import yaml

THRESH = 99.9  # per-build threshold (gate-threshold branch); 100 before
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
SETS = {
    "LCG": f"{R}/2026-09-28_lcg_holdout/runs2/*/detection_report.yaml",
    "Jena": f"{R}/2026-09-28_mucor_jena_holdout/runs_scaffolds/*/detection_report.yaml",
    "Mucoro293": f"{R}/2026-09-28_next_fixes/Mucoromycota_a9fb69c/runs/*/detection_report.yaml",
    "Zygo23_scaffold": f"{R}/2026-09-28_next_fixes/zygo23_a9fb69c/scaffold/runs/*/detection_report.yaml",
    "Zygo23_contig": f"{R}/2026-09-28_next_fixes/zygo23_a9fb69c/contig/runs/*/detection_report.yaml",
}


def modelled(g):
    # as flank_carried._modelled: polished by a tool, or annotated (proteome fast path)
    return str(g.get("status", "")).startswith("polished") or g.get("method") == "diamond_proteome"


def aa_len(g):
    if g.get("exons"):
        return sum(e["end"] - e["start"] + 1 for e in g["exons"]) // 3
    return (abs(g["end"] - g["start"]) + 1) // 3  # fast-path: span as an upper bound


rows = []
for sname, pat in SETS.items():
    for f in sorted(glob.glob(pat)):
        org = f.split("/")[-2]
        try:
            r = yaml.safe_load(open(f))
        except Exception:
            continue
        for ci, c in enumerate(r.get("detected") or []):
            if not c["family"].startswith("Mucoromycota"):
                continue
            clf = c.get("idiomorph_classifier") or {}
            idio = c["idiomorph"]
            scores = clf.get("scores") or {}
            top = max(scores.values()) if scores else None
            inp = clf.get("classifier_input")
            ge = c.get("gene_evidence") or []
            flanks = sorted({g["gene"] for g in ge if g["role"].startswith("flanking") and modelled(g)
                             and g.get("identity") is not None and g["identity"] >= 40})
            want = {"Minus": "sexM", "Plus": "sexP"}.get(idio)
            core = [g for g in ge if g["gene"] == want and modelled(g)]
            clen = max((aa_len(g) for g in core), default=0)
            if inp != "model":
                clen = 0
            if not clf:
                route = "no_classifier"
            elif c.get("split_locus"):
                route = "split_locus"
            elif inp == "model" and top is not None and top >= THRESH:
                route = "score"
            else:
                route = "flank"
            r1 = "rnhA" in flanks
            r2a, r2b = clen >= 100, clen >= 125
            rows.append(dict(set=sname, org=org, call=ci, idiomorph=idio, confidence=c["confidence"],
                             locus_class=c.get("locus_class"), clf_input=inp, top_score=top,
                             margin=clf.get("margin"), core_len=clen, flanks=",".join(flanks),
                             route=route,
                             drop_R1=route == "flank" and not r1,
                             drop_R2a=route == "flank" and not r2a,
                             drop_R2b=route == "flank" and not r2b,
                             drop_R3=route == "flank" and not (r1 or r2a)))
with open("replay_rules.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
for s in SETS:
    rs = [x for x in rows if x["set"] == s]
    genomes = {x["org"] for x in rs}
    line = f"{s}: calls {len(rs)} genomes {len(genomes)} flank-route {sum(x['route']=='flank' for x in rs)}"
    for k in ("R1", "R2a", "R2b", "R3"):
        d = [x for x in rs if x[f"drop_{k}"]]
        left = {x["org"] for x in rs if not x[f"drop_{k}"]}
        line += f" | {k} drop {len(d)} (genomes losing all calls {len(genomes - left)})"
    print(line)
