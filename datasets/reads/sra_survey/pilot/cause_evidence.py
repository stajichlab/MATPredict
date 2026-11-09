#!/usr/bin/env python3
"""Evidence for why a reads species cannot be typed by presence, or why its panel failed (#195, #196): one row per species, joining what the
assemblies say about its MAT loci with what the reads gate and the panel build said. The curator's decisions are kept separately in causes.tsv
(species, cause, basis), which wave_driver.py reads; this table is the evidence those decisions cite.

    cause_evidence.py --pilot PILOT_DIR --release RESULTS_DIR [--release RESULTS_DIR ...] [--plans WAVES_DIR] --out cause_evidence.tsv [--species NAME ...]

PILOT_DIR: committed wave folders (wave*/gates_final.tsv). RESULTS_DIR: a MATPredict release folder with genomes.tsv and loci.tsv (for example
results/2026-11-05_v081_ascomycota). WAVES_DIR (HPCC): the wave folders with plan/<slug>/{status, panel/*_lineages.tsv}; without it the build columns
are blank. Default species: every species whose decision is PROVISIONAL_HOMOTHALLIC, a ZYGOSITY word, released_*, NOT_TYPABLE, or a HOLD from the
panel stage (panel_build_failed, panel_L1_error, panel_no_unique_kmers).

Assembly columns count genomes of the species (by genomes.tsv `species`):
  genomes, genomes_with_mat (a called locus), single (one idiomorph label), both_one_locus (a homothallic_candidate locus: both genes in one locus),
  both_separate (two or more idiomorph labels at separate loci: silent cassettes, unlinked homothallic loci, or a collapsed diploid), labels (label: genomes).
Build columns: plan_status, lineages ("<loci> loci: <a> <name1> / <b> <name2>" per lineage, from the panel's lineage table).
Standard library only.
"""
import argparse, csv, glob, os, re
from collections import Counter, defaultdict

BUILD_HOLDS = ("panel_build_failed", "panel_L1_error", "panel_no_unique_kmers")


def rd(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def gates(pilot):
    """{slug: (species, decision, reasons)} from every committed wave (the later wave wins for a slug that appears twice)."""
    g = {}
    for w in sorted(glob.glob(os.path.join(pilot, "wave*"))):
        f = os.path.join(w, "gates_final.tsv")
        if re.fullmatch(r"wave\d+[a-z_]*", os.path.basename(w)) and os.path.exists(f):
            for r in rd(f): g[r["slug"]] = (r["species"], r["decision"], r["reasons"])
    return g


def is_target(decision):
    return (decision.startswith(("PROVISIONAL_HOMOTHALLIC", "PASS_ZYGOSITY", "PROVISIONAL_ZYGOSITY", "NOT_TYPABLE")) or "zygosity" in decision
            or "released_" in decision or any(h in decision for h in BUILD_HOLDS))


def assembly_evidence(release_dirs):
    """{species: Counter of the assembly columns plus a Counter of labels} from the release tables."""
    sp_of = {}
    for d in release_dirs:
        for r in rd(os.path.join(d, "genomes.tsv")): sp_of[r["genome"]] = r["species"]
    per = defaultdict(lambda: {"labels": set(), "homo": False})
    for d in release_dirs:
        for r in rd(os.path.join(d, "loci.tsv")):
            g = per[r["genome"]]
            if r["locus_class"] == "homothallic_candidate": g["homo"] = True
            if r["idiomorph"] not in ("", "undetermined"): g["labels"].add(r["idiomorph"])
    out = defaultdict(lambda: {"n": Counter(), "labels": Counter()})
    for genome, species in sp_of.items():
        e = out[species]; e["n"]["genomes"] += 1
        g = per.get(genome)
        if not g or not (g["labels"] or g["homo"]): continue
        e["n"]["genomes_with_mat"] += 1
        if g["homo"]: e["n"]["both_one_locus"] += 1
        elif len(g["labels"]) >= 2: e["n"]["both_separate"] += 1
        else: e["n"]["single"] += 1
        for lab in g["labels"]: e["labels"][lab] += 1
    return out


def build_evidence(waves, slug):
    """(plan status, lineage summary) of a species' panel build, from the last wave that planned it."""
    ds = sorted(glob.glob(os.path.join(waves, "*", "plan", slug)))
    if not ds: return "", ""
    d = ds[-1]
    st = open(os.path.join(d, "status")).read().strip() if os.path.exists(os.path.join(d, "status")) else ""
    lt = glob.glob(os.path.join(d, "panel", "*_lineages.tsv"))
    if not lt: return st, "no lineage table"
    rows = rd(lt[0])
    names = [c for c in rows[0].keys() if c not in ("lineage", "loci", "kept") and not c.startswith("rep_")] if rows else []
    lin = "; ".join(f"L{r['lineage']} {r['loci']} loci: " + " / ".join(f"{r[n]} {n}" for n in names) + ("" if r.get("kept") in ("yes", "") else " (dropped)") for r in rows)
    return st, lin or "no lineages"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pilot", required=True); ap.add_argument("--release", action="append", required=True); ap.add_argument("--plans")
    ap.add_argument("--out", required=True); ap.add_argument("--species", nargs="*")
    a = ap.parse_args()
    g = gates(a.pilot); ev = assembly_evidence(a.release)
    cols = ["slug", "species", "decision", "reads", "genomes", "genomes_with_mat", "single", "both_one_locus", "both_separate", "labels", "plan_status", "lineages"]
    rows = []
    for slug, (species, decision, reasons) in sorted(g.items(), key=lambda kv: kv[1][0]):
        if a.species: keep = species in a.species
        else: keep = is_target(decision)
        if not keep: continue
        e = ev.get(species, {"n": Counter(), "labels": Counter()})
        m = re.search(r"done=(\d+) typed=(\d+) both=(\d+)", reasons)
        st, lin = build_evidence(a.plans, slug) if a.plans else ("", "")
        rows.append({"slug": slug, "species": species, "decision": decision, "reads": f"done {m[1]}, single {m[2]}, both {m[3]}" if m else "",
                     **{k: e["n"][k] for k in ("genomes", "genomes_with_mat", "single", "both_one_locus", "both_separate")},
                     "labels": ", ".join(f"{k}: {v}" for k, v in e["labels"].most_common()), "plan_status": st, "lineages": lin})
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"{len(rows)} species -> {a.out}")


if __name__ == "__main__":
    main()
