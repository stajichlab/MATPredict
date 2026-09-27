"""Before/after for the Russulaceae PR record on the 15 Russulales genomes.

before = results/2026-09-27_receptor_curation/run_29818af (Heterobasidion PR record only)
after  = run_4a0f0eb (+ Russula nobilis 2830151_kdtol00553_PR_B1)
Same code; only the database differs. Reports PR/HD calls per genome and genus,
whether PR calls overlap a strict-CAAX flagged STE3 locus, PR vs HD contig,
withheld PR loci, and wall time.
"""
import csv
import os
import statistics
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
BEFORE = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_receptor_curation/run_29818af"
AFTER = f"{HERE}/run_4a0f0eb"
POS = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_pheromone_positional/uncurated_loci.tsv"
samples = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
flag = {}
for r in csv.DictReader(open(POS), delimiter="\t"):
    if r["flag_strict"] == "1":
        flag.setdefault(r["asm"], []).append((r["contig"], int(r["start"]), int(r["end"])))
panel = [l.split("\t")[0] for l in open(f"{HERE}/lists/russulales.tsv")]


def load(root, a):
    p = f"{root}/runs/{a}/detection_report.yaml"
    if not os.path.exists(p):
        return None
    r = yaml.safe_load(open(p)) or {}
    w = f"{root}/runs/{a}/wall_seconds"
    r["_wall"] = float(open(w).read().strip()) if os.path.exists(w) else None
    return r


def calls(r, fam):
    return [x for x in (r.get("detected") or []) if x["family"] == fam] if r else []


def withheld(r, fam):
    return [x for x in (r.get("suppressed_loci") or []) if x["family"] == fam] if r else []


rows, wb, wa = [], [], []
for a in panel:
    s = samples.get(a, {})
    rb, ra = load(BEFORE, a), load(AFTER, a)
    prb, pra = calls(rb, "Basidiomycota:PR"), calls(ra, "Basidiomycota:PR")
    hda = calls(ra, "Basidiomycota:HD")
    fl = flag.get(a, [])
    at_flag = sum(any(c == x["contig"] and s0 <= x["end"] and e0 >= x["start"] for c, s0, e0 in fl) for x in pra)
    shared = sum(x["contig"] in {h["contig"] for h in hda} for x in pra)
    ev = "; ".join(f"{e['gene']}:{e.get('identity')}:{e.get('status')}"
                   for x in pra for e in x.get("gene_evidence", []))
    if rb and rb.get("_wall"):
        wb.append(rb["_wall"])
    if ra and ra.get("_wall"):
        wa.append(ra["_wall"])
    rows.append(dict(asm=a, family=s.get("FAMILY"), genus=s.get("GENUS"), species=s.get("SPECIES"),
                     PR_before=len(prb), PR_after=len(pra), HD_after=len(hda),
                     PR_at_flagged=at_flag, PR_on_HD_contig=shared,
                     PR_withheld_after=len(withheld(ra, "Basidiomycota:PR")),
                     PR_calls=";".join(f"{x['contig']}:{x['start']}-{x['end']}/{x['confidence']}/{x['locus_class']}"
                                       for x in pra),
                     PR_evidence=ev, wall_before=rb.get("_wall") if rb else None,
                     wall_after=ra.get("_wall") if ra else None))
with open(f"{HERE}/per_genome.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"genomes {len(rows)}")
print(f"PR genomes called before {sum(r['PR_before'] > 0 for r in rows)}, after {sum(r['PR_after'] > 0 for r in rows)}")
print(f"PR calls after {sum(r['PR_after'] for r in rows)}; at flagged STE3 locus {sum(r['PR_at_flagged'] for r in rows)}; "
      f"on an HD contig {sum(r['PR_on_HD_contig'] for r in rows)}")
print(f"HD genomes called after {sum(r['HD_after'] > 0 for r in rows)}")
print(f"median wall before {statistics.median(wb):.0f} s, after {statistics.median(wa):.0f} s")
by = {}
for r in rows:
    k = f"{r['family']}/{r['genus']}"
    b, c = by.get(k, (0, 0))
    by[k] = (b + 1, c + (r["PR_after"] > 0))
for k, (n, c) in sorted(by.items()):
    print(f"  {k:40s} {c}/{n} PR called")
for r in rows:
    if r["PR_after"] or r["PR_before"]:
        print(f"  CALL {r['species']:35s} {r['PR_before']}->{r['PR_after']} {r['PR_calls']} | {r['PR_evidence']}")
