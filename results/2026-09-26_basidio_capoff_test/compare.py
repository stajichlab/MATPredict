"""Compare cap-off re-runs with the capped full run, same code (run-ad1f865).

For each selected genome: calls with cap off vs capped, wall time both ways,
and for every gained call its gene evidence (identity, coverage, status) and
whether it overlaps a locus the capped run withheld (suppressed_loci).
Writes gained_calls.tsv, per_genome.tsv; prints a summary.
"""
import csv, os, statistics, yaml

O = os.path.dirname(os.path.abspath(__file__))
B = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_basidiomycota_full"
sel = list(csv.DictReader(open(f"{O}/selected.tsv"), delimiter="\t"))


def load(d):
    f = f"{d}/detection_report.yaml"
    if not os.path.exists(f):
        return None, None
    w = f"{d}/wall_seconds"
    wall = float(open(w).read().strip()) if os.path.exists(w) else None
    return yaml.safe_load(open(f)), wall


def find_off(g):
    for ch in sorted(os.listdir(O)):
        d = f"{O}/{ch}/runs/{g}"
        if ch.startswith("capoff_") and os.path.isdir(d):
            return d
    return None


def overlap(a, b):
    return a["contig"] == b["contig"] and a["start"] <= b["end"] and b["start"] <= a["end"]


pg, gained, walls = [], [], []
for s in sel:
    g = s["genome"]
    on, won = load(f"{B}/{s['chunk_capped']}/runs/{g}")
    doff = find_off(g)
    off, woff = load(doff) if doff else (None, None)
    con = (on or {}).get("detected") or []
    coff = (off or {}).get("detected") or [] if off else []
    status = "no_report" if off is None else ("gained" if len(coff) > len(con) else "same")
    pg.append(dict(genome=g, order=s["order"], species=s["species"], capped_calls=len(con),
                   capoff_calls=len(coff) if off else "", wall_capped=won, wall_capoff=woff,
                   status=status))
    if won and woff:
        walls.append((won, woff))
    if not off:
        continue
    for c in coff:
        if any(overlap(c, x) for x in con):
            continue
        sup = [x for x in (on.get("suppressed_loci") or []) if overlap(c, x)]
        ev = "; ".join(f"{e['gene']}:{e['identity']}/{e.get('coverage')}/{e['status']}"
                       for e in c.get("gene_evidence", []))
        core = [e for e in c.get("gene_evidence", []) if e["role"] == "core_MAT"
                and str(e["status"]).startswith("polished")]
        gained.append(dict(genome=g, order=s["order"], species=s["species"], family=c["family"],
                           contig=c["contig"], start=c["start"], end=c["end"],
                           idiomorph=c["idiomorph"], confidence=c["confidence"],
                           locus_class=c["locus_class"], polished_genes=c.get("polished_genes"),
                           modelled_core=len(core),
                           best_core_identity=max([e["identity"] for e in core], default=""),
                           overlaps_capped_withheld=bool(sup), evidence=ev))

for name, rows in (("per_genome.tsv", pg), ("gained_calls.tsv", gained)):
    if rows:
        with open(f"{O}/{name}", "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)

n_rep = sum(1 for r in pg if r["status"] != "no_report")
print(f"genomes selected {len(pg)}; cap-off reports {n_rep}; genomes gaining a call "
      f"{sum(1 for r in pg if r['status'] == 'gained')}; gained calls {len(gained)}")
if walls:
    ratio = [b / a for a, b in walls if a]
    print(f"wall capped median {statistics.median(a for a, _ in walls):.0f} s, cap-off median "
          f"{statistics.median(b for _, b in walls):.0f} s; ratio median {statistics.median(ratio):.2f}, "
          f"max {max(ratio):.2f}; total capped {sum(a for a, _ in walls)/3600:.1f} h, "
          f"cap-off {sum(b for _, b in walls)/3600:.1f} h")
for r in gained:
    print("GAINED", r["order"], r["species"], r["family"], r["idiomorph"], r["confidence"],
          r["locus_class"], "modelled_core", r["modelled_core"], "best_core_id", r["best_core_identity"],
          "overlaps_withheld", r["overlaps_capped_withheld"])
    print("    ", r["evidence"])
