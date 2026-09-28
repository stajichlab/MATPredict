"""Review-fix measurement (branch review-fixes, ca3045c).

Mucoromycota: new results/2026-09-28_review_fixes/Mucoromycota_f25cf70 vs
              results/2026-09-27_split_locus/Mucoromycota_4174440 (8d80bed+split, PR #9 44567fc code).
Agaricales:   new results/2026-09-28_review_fixes/agaricales_panel
              (ca3045c + basidio-anchors db) vs results/2026-09-27_locus_merge/agaricales_panel.
"""
import collections, csv, glob, yaml

BASE = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
samp = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}


def load(d):
    out = {}
    for f in glob.glob(f"{d}/runs/*/detection_report.yaml"):
        try:
            out[f.split("/")[-2]] = yaml.safe_load(open(f)) or {}
        except Exception:
            pass
    return out


def calls(rep):
    return [(x["family"], x["contig"], x["start"], x["end"], x["idiomorph"], x["confidence"])
            for x in rep.get("detected") or []]


def same(a, b):
    return a[0] == b[0] and a[1] == b[1] and a[2] <= b[3] and a[3] >= b[2]


# ---------------- Mucoromycota
old = load(f"{BASE}/2026-09-27_split_locus/Mucoromycota_4174440")
new = load(f"{BASE}/2026-09-28_review_fixes/Mucoromycota_f25cf70")
both = sorted(set(old) & set(new))
print(f"== Mucoromycota: genomes old {len(old)} new {len(new)} both {len(both)}")
gained, lost, changed = [], [], []
floor_rows = collections.Counter()
floor_ex = []
for g in both:
    o, n = calls(old[g]), calls(new[g])
    for c in n:
        m = [x for x in o if same(x, c)]
        if not m:
            gained.append((g, c))
        elif (m[0][4], m[0][5]) != (c[4], c[5]):
            changed.append((g, m[0][4], m[0][5], c[4], c[5], c[0]))
    for c in o:
        if not any(same(c, x) for x in n):
            lost.append((g, c))
    for s in new[g].get("suppressed_loci") or []:
        floor_rows[s.get("withheld_reason")] += 1
        if s.get("withheld_reason") == "below_fraction_floor" and len(floor_ex) < 8:
            floor_ex.append((g, samp.get(g, {}).get("SPECIES"), s["contig"], s["start"],
                             s.get("fraction_found"), s.get("best_identity")))
called_old = sum(1 for g in both if old[g].get("detected"))
called_new = sum(1 for g in both if new[g].get("detected"))
print(f" genomes called {called_old} -> {called_new}; calls gained {len(gained)} lost {len(lost)} changed {len(changed)}")
for g, c in gained:
    sp = samp.get(g, {}).get("SPECIES")
    rp = [x for x in new[g]["detected"] if x["contig"] == c[1] and x["start"] == c[2]][0]
    print(f"   GAINED {sp} {g} {c[0]} {c[1]}:{c[2]} {c[4]}/{c[5]} pass={rp.get('detection_pass')} polished={rp.get('polished_genes')}")
for g, c in lost:
    print(f"   LOST {samp.get(g, {}).get('SPECIES')} {g} {c}")
for x in changed:
    print(f"   CHANGED {samp.get(x[0], {}).get('SPECIES')} {x}")
print(" suppressed_loci by reason:", dict(floor_rows))
for x in floor_ex:
    print("   floor:", x)
cls = sum(1 for g in both for s in new[g].get("suppressed_loci") or [] if s.get("idiomorph_classifier"))
print(" suppressed loci carrying a classifier block:", cls)
for g in both:
    sp = samp.get(g, {}).get("SPECIES", "")
    if g.startswith(("GCA_900175165", "GCA_900079185")):
        print(f"   F1 case {sp} {g}: {calls(new[g])}")

# ---------------- Agaricales
old = load(f"{BASE}/2026-09-27_locus_merge/agaricales_panel")
new = load(f"{BASE}/2026-09-28_review_fixes/agaricales_panel")
both = sorted(set(old) & set(new))
print(f"\n== Agaricales panel: genomes old {len(old)} new {len(new)} both {len(both)}")
fo, fn = collections.Counter(), collections.Counter()
merges_o, merges_n = collections.Counter(), collections.Counter()
sub_present, sub_labels, conf_changes, unv_o, unv_n = 0, collections.Counter(), [], 0, 0
for g in both:
    for x in old[g].get("detected") or []:
        fo[x["family"]] += 1
        if x.get("merged_from"):
            merges_o["+".join(sorted(m["family"].split(":")[1] for m in x["merged_from"]))] += 1
        unv_o += (x.get("verification") or {}).get("status") == "unverified"
    for x in new[g].get("detected") or []:
        fn[x["family"]] += 1
        unv_n += (x.get("verification") or {}).get("status") == "unverified"
        if x.get("merged_from"):
            merges_n["+".join(sorted(m["family"].split(":")[1] for m in x["merged_from"]))] += 1
            if x.get("subloci"):
                sub_present += 1
                for s in x["subloci"]:
                    sub_labels[s["sublocus"]] += 1
    oc = {(x["family"], x["contig"], x["start"]): x["confidence"] for x in old[g].get("detected") or []}
    for x in new[g].get("detected") or []:
        k = (x["family"], x["contig"], x["start"])
        if k in oc and oc[k] != x["confidence"]:
            conf_changes.append((g, k, oc[k], x["confidence"]))
print(" calls", sum(fo.values()), "->", sum(fn.values()))
for k in sorted(set(fo) | set(fn)):
    print(f"   {k}: {fo[k]} -> {fn[k]}")
print(" merges old", dict(merges_o))
print(" merges new", dict(merges_n))
print(" merged calls carrying subloci:", sub_present, "of", sum(merges_n.values()), dict(sub_labels))
print(" unverified calls", unv_o, "->", unv_n)
print(" confidence changes on same call:", len(conf_changes))
for c in conf_changes[:10]:
    print("   ", c)
