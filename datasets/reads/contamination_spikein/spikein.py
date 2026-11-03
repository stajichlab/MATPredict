#!/usr/bin/env python3
"""Contamination spike-in for the raw-reads screen: how little of a second organism can GmG see in a read set?

Real reads of a background organism (for example a clean *Coccidioides*) are mixed with real reads of a
contaminant (*Histoplasma*, *Lecanicillium*, or the other *Coccidioides* species) at known fractions. GmG then
screens the first N pairs of each mixed set at several N, so one file answers "is 2M pairs enough, or are 8M needed".

  spikein.py mix    --bg R1 R2 --con R1 R2 --fraction F --pairs N --seed S --out-r1 X --out-r2 Y
                    Write N pairs; each pair is the contaminant's with probability F (seeded), else the background's.
                    The pairs of each source are used in file order, the mix is random, so ANY prefix of the output is
                    a mix at about F. A sidecar X.mix.json records the realised contaminant pairs at each checkpoint.
  spikein.py sheet  --sources TSV --outdir DIR [--fractions ...]
                    Write conditions.tsv and the GmG sample sheet for the mixed files (one id, two rows R1 and R2).
  spikein.py mix-one --sources TSV --outdir DIR --index I [--pairs N]
                    Build the mixed files of condition I (1-based row of conditions.tsv); a SLURM array runs this.
                    Skipped when its sidecar already exists, so a requeue is safe.
  spikein.py summarize --conditions TSV --runs DEPTH=results.jsonl ... [--clean-jsonl results.jsonl] --out PREFIX
                    One row per condition and depth with GmG's call and its mixture block, and a detection table.

GmG counts --max-reads per file, so --max-reads N on a sheet of paired files is N pairs.
"""
import argparse
import csv
import gzip
import json
import os
import random
import sys

CHECKPOINTS = (100_000, 500_000, 2_000_000, 6_000_000)
FRACTIONS = (0.005, 0.01, 0.02, 0.05, 0.10, 0.25)


# ---------------------------------------------------------------------------------------------- mixing
def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_pairs(r1, r2):
    """Yield (record1, record2) FASTQ records (four lines each, as text) from two paired files in file order."""
    with _open(r1) as f1, _open(r2) as f2:
        while True:
            a = [f1.readline() for _ in range(4)]
            b = [f2.readline() for _ in range(4)]
            if not a[0] or not b[0]:
                return
            yield "".join(a), "".join(b)


def mix(bg, con, fraction, n_pairs, seed, out_r1, out_r2, checkpoints=CHECKPOINTS):
    """Write `n_pairs` pairs to out_r1 / out_r2. Pair i comes from the contaminant `con` with probability `fraction`
    (a seeded draw), otherwise from the background `bg`; each source is read in file order. bg and con are (r1, r2).
    Returns {"n_pairs", "fraction", "contaminant_pairs", "checkpoints": {N: contaminant pairs among the first N}}."""
    rng = random.Random(seed)
    bg_it, con_it = read_pairs(*bg), read_pairs(*con)
    marks = sorted(c for c in checkpoints if c <= n_pairs)
    seen = {}
    from_con = 0
    with _open(out_r1, "wt") as o1, _open(out_r2, "wt") as o2:
        for i in range(1, n_pairs + 1):
            use_con = fraction > 0 and rng.random() < fraction
            try:
                a, b = next(con_it if use_con else bg_it)
            except StopIteration:
                raise SystemExit(f"the {'contaminant' if use_con else 'background'} reads ran out after {i - 1} mixed pairs "
                                 f"(wanted {n_pairs} at fraction {fraction})")
            o1.write(a)
            o2.write(b)
            from_con += use_con
            if i in marks:
                seen[i] = from_con
    return {"n_pairs": n_pairs, "fraction": fraction, "seed": seed, "contaminant_pairs": from_con, "checkpoints": seen}


# ---------------------------------------------------------------------------------------------- sheet
def read_sources(path):
    """sources.tsv: id, r1, r2, species, role (background or contaminant), and optional group and only_with
    (a contaminant used with only these backgrounds, comma separated)."""
    return list(csv.DictReader(open(path), delimiter="\t"))


def conditions(sources, fractions=FRACTIONS):
    """Every background x contaminant x fraction, plus a fraction-0 control per background. A same-genus pair
    (group equal) is kept: it is the negative that a genus-level screen should not flag."""
    bgs = [s for s in sources if s["role"] == "background"]
    cons = [s for s in sources if s["role"] == "contaminant"]
    out = []
    for b in bgs:
        out.append({"bg": b["id"], "con": "none", "fraction": 0.0, "kind": "control"})
        for c in cons:
            only = [x for x in (c.get("only_with") or "").split(",") if x]
            if only and b["id"] not in only:
                continue
            kind = "same_genus" if b.get("group") and b.get("group") == c.get("group") else "other_genus"
            for f in fractions:
                out.append({"bg": b["id"], "con": c["id"], "fraction": f, "kind": kind})
    for i, c in enumerate(out):
        c["id"] = f"{c['bg']}__{c['con']}__f{c['fraction']:g}"
        c["seed"] = 1000 + i
    return out


# ---------------------------------------------------------------------------------------------- summary
def gmg_row(r):
    """The fields of one GmG results.jsonl record that matter here."""
    q = r.get("qc") or {}
    m = q.get("mixture") or {}
    s = m.get("second") or {}
    calls = (r.get("guess") or {}).get("calls") or {}
    genus = (calls.get("genus") or {}).get("taxon")
    species = (calls.get("species") or {}).get("taxon") or (calls.get("species") or {}).get("vote")
    return {
        "id": r.get("id"), "status": r.get("status"), "best_rank": r.get("best_rank"), "genus": genus, "species": species,
        "flags": ",".join(r.get("qc_flags") or []), "mixture_suspected": bool(m.get("suspected")),
        "coverage_est": q.get("coverage_est"), "n_hashes": q.get("n_hashes"), "primary": m.get("primary"),
        "explained": m.get("explained_fraction"), "second": s.get("taxon"), "second_hashes": s.get("hashes"),
        "second_explained": s.get("explained_fraction"), "second_ref_coverage": s.get("ref_coverage"),
        "second_concentration": s.get("concentration"),
    }


def read_jsonl(path):
    rows = []
    for line in open(path):
        line = line.strip()
        if line:
            rows.append(gmg_row(json.loads(line)))
    return rows


def sees(row, contaminant_genus, min_hashes, min_concentration):
    """True when the row names the contaminant genus as the second organism with enough hashes and enough of them
    concentrated in one reference genome (a real contaminant's hashes sit in one genome; noise hashes are scattered),
    or as the primary organism (the contaminant dominates the set)."""
    if row["primary"] == contaminant_genus and row["second"] and row["second"] != contaminant_genus:
        return True
    return bool(row["second"] == contaminant_genus and (row["second_hashes"] or 0) >= min_hashes
                and (row["second_concentration"] or 0) >= min_concentration)


def sensitivity_table(rows, conds, genus_of, thresholds):
    """Per (contaminant, fraction, depth): how many of the conditions are seen by GmG's own flag and by each
    (min_hashes, min_concentration) rule. `rows` are summary rows with cond fields; genus_of maps a contaminant id
    to its genus. Returns a list of dicts."""
    out = []
    keys = sorted({(r["con"], r["fraction"], r["depth"]) for r in rows if r["kind"] == "other_genus"})
    for con, f, depth in keys:
        sel = [r for r in rows if (r["con"], r["fraction"], r["depth"]) == (con, f, depth)]
        g = genus_of[con]
        rec = {"contaminant": con, "fraction": f, "depth_pairs": depth, "n": len(sel),
               "gmg_flag": sum(1 for r in sel if r["mixture_suspected"] and (r["second"] == g or r["primary"] == g))}
        for h, c in thresholds:
            rec[f"seen_h{h}_q{c:g}"] = sum(1 for r in sel if sees(r, g, h, c))
        out.append(rec)
    return out


def false_positives(neg_rows, thresholds):
    """Rows that should show nothing (fraction-0 controls and clean sets): how many a rule would flag, whatever the
    second genus is. Returns {"n", "gmg_flag", "h<h>_c<c>": count}."""
    out = {"n": len(neg_rows), "gmg_flag": sum(1 for r in neg_rows if r["mixture_suspected"])}
    for h, c in thresholds:
        out[f"h{h}_q{c:g}"] = sum(1 for r in neg_rows if r["second"] and (r["second_hashes"] or 0) >= h and (r["second_concentration"] or 0) >= c)
    return out


# ---------------------------------------------------------------------------------------------- CLI
def _cmd_mix(a):
    meta = mix(tuple(a.bg), tuple(a.con), a.fraction, a.pairs, a.seed, a.out_r1, a.out_r2)
    with open(a.out_r1 + ".mix.json", "w") as fh:
        json.dump(meta, fh)
    print(json.dumps(meta))


def _cmd_mixone(a):
    sources = {s["id"]: s for s in read_sources(a.sources)}
    conds = list(csv.DictReader(open(os.path.join(a.outdir, "conditions.tsv")), delimiter="\t"))
    c = conds[a.index - 1]
    mixed = os.path.join(a.outdir, "mixed")
    os.makedirs(mixed, exist_ok=True)
    o1, o2 = os.path.join(mixed, c["id"] + "_R1.fastq.gz"), os.path.join(mixed, c["id"] + "_R2.fastq.gz")
    side = os.path.join(mixed, c["id"] + ".mix.json")
    if os.path.exists(side):
        print("already built:", c["id"])
        return
    bg = sources[c["bg"]]
    con = sources[c["con"]] if c["con"] != "none" else bg      # a control never draws from its second source
    p1, p2 = o1.replace(".fastq.gz", ".part.fastq.gz"), o2.replace(".fastq.gz", ".part.fastq.gz")      # still ends in .gz: that is how mix() chooses gzip
    meta = mix((bg["r1"], bg["r2"]), (con["r1"], con["r2"]), float(c["fraction"]), a.pairs, int(c["seed"]), p1, p2)
    os.replace(p1, o1)
    os.replace(p2, o2)
    with open(side, "w") as fh:        # the sidecar is written last: its presence means the files are complete
        json.dump({**meta, "id": c["id"], "bg": c["bg"], "con": c["con"]}, fh)
    print(json.dumps(meta))


def _cmd_sheet(a):
    sources = read_sources(a.sources)
    conds = conditions(sources, a.fractions)
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "conditions.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "bg", "con", "fraction", "kind", "seed"], delimiter="\t")
        w.writeheader()
        w.writerows(conds)
    byid = {s["id"]: s for s in sources}
    with open(os.path.join(a.outdir, "samples.tsv"), "w") as fh:
        fh.write("id\tpath\texpected\n")
        for c in conds:
            for n in (1, 2):
                fh.write(f"{c['id']}\t{os.path.join(a.outdir, 'mixed', c['id'] + f'_R{n}.fastq.gz')}\t{byid[c['bg']]['species']}\n")
    print(f"{len(conds)} conditions -> {a.outdir}")


def _cmd_summarize(a):
    conds = {c["id"]: c for c in csv.DictReader(open(a.conditions), delimiter="\t")}
    sources = {s["id"]: s for s in read_sources(a.sources)}
    genus_of = {k: v["species"].split()[0] for k, v in sources.items()}
    thresholds = [(h, q) for h in (3, 5, 10, 20, 50) for q in (0.5, 0.7, 0.8)]
    rows, neg = [], []
    for spec in a.runs:
        depth, path = spec.split("=", 1)
        for r in read_jsonl(path):
            cond = conds.get(r["id"])
            if not cond:
                continue
            row = {**r, "depth": int(depth), "bg": cond["bg"], "con": cond["con"], "fraction": float(cond["fraction"]), "kind": cond["kind"]}
            rows.append(row)
            if cond["kind"] == "control":
                neg.append(row)
    clean = read_jsonl(a.clean_jsonl) if a.clean_jsonl else []
    cols = ["id", "depth", "bg", "con", "fraction", "kind", "status", "best_rank", "genus", "species", "flags", "mixture_suspected",
            "coverage_est", "n_hashes", "primary", "explained", "second", "second_hashes", "second_explained", "second_ref_coverage", "second_concentration"]
    with open(a.out + ".rows.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["bg"], r["con"], r["fraction"], r["depth"])))
    sens = sensitivity_table(rows, conds, genus_of, thresholds)
    with open(a.out + ".detection.tsv", "w") as fh:
        if sens:
            w = csv.DictWriter(fh, fieldnames=list(sens[0]), delimiter="\t")
            w.writeheader()
            w.writerows(sens)
    same = [r for r in rows if r["kind"] == "same_genus"]
    same_summary = {"n": len(same), "gmg_flag": sum(1 for r in same if r["mixture_suspected"]),
                    "species_call_changed": sum(1 for r in same if r["species"] and r["species"] != sources[r["bg"]]["species"]),
                    "any_other_genus_second": sum(1 for r in same if r["second"] and r["second"] != genus_of[r["bg"]])}
    fp = {"controls": false_positives(neg, thresholds), "clean_sets_2M": false_positives(clean, thresholds) if clean else None,
          "same_genus_mixtures": same_summary}
    json.dump(fp, open(a.out + ".false_positives.json", "w"), indent=1)
    print(f"{len(rows)} rows, {len(sens)} detection cells; false positives: {json.dumps(fp)}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("mix")
    m.add_argument("--bg", nargs=2, required=True)
    m.add_argument("--con", nargs=2, required=True)
    m.add_argument("--fraction", type=float, required=True)
    m.add_argument("--pairs", type=int, default=8_000_000)
    m.add_argument("--seed", type=int, default=1)
    m.add_argument("--out-r1", required=True)
    m.add_argument("--out-r2", required=True)
    m.set_defaults(func=_cmd_mix)
    o = sub.add_parser("mix-one")
    o.add_argument("--sources", required=True)
    o.add_argument("--outdir", required=True)
    o.add_argument("--index", type=int, required=True)
    o.add_argument("--pairs", type=int, default=6_000_000)
    o.set_defaults(func=_cmd_mixone)
    s = sub.add_parser("sheet")
    s.add_argument("--sources", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--fractions", type=float, nargs="+", default=list(FRACTIONS))
    s.set_defaults(func=_cmd_sheet)
    z = sub.add_parser("summarize")
    z.add_argument("--conditions", required=True)
    z.add_argument("--sources", required=True)
    z.add_argument("--runs", nargs="+", required=True, help="DEPTH=results.jsonl, one per GmG run (DEPTH in pairs)")
    z.add_argument("--clean-jsonl")
    z.add_argument("--out", required=True)
    z.set_defaults(func=_cmd_summarize)
    a = p.parse_args(argv)
    a.func(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
