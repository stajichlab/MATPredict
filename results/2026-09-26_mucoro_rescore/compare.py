"""Label changes between the old Mucoromycota scan and a re-score.

usage: compare.py OLD_RUNS NEW_RUNS [--sets]

Calls are paired within a genome by contig overlap. Reports, overall and per
family: unchanged, Plus->Minus, Minus->Plus, other label changes, new calls,
lost calls. With --sets, also the three watched groups:
  * the 33 Rhizopus sexM-only calls labelled Plus in OLD
  * the 13 Minus calls whose HMG box sits in the sexP clade
    (results/2026-09-26_sexMP_phylogeny/candidates.tsv)
  * every Lichtheimiaceae / Syncephalastraceae genome
and the relaxed-pass calls (evidence diagnostics `relaxed_call` rows).
"""
import collections, csv, glob, json, os, sys
import yaml
from yaml import CSafeLoader as L

OLD, NEW = sys.argv[1], sys.argv[2]
SETS = "--sets" in sys.argv
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
CAND = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_sexMP_phylogeny/candidates.tsv"
meta = {r["ASMID"]: r for r in csv.DictReader(open(SAMPLES))}


def calls(runs, g):
    p = f"{runs}/{g}/detection_report.yaml"
    if not os.path.exists(p):
        return None
    return (yaml.load(open(p), Loader=L) or {}).get("detected") or []


def ov(a, b):
    return a["contig"] == b["contig"] and not (a["end"] < b["start"] or a["start"] > b["end"])


genomes = sorted(set(os.listdir(OLD)) & set(os.listdir(NEW)))
tot = collections.Counter()
fam = collections.defaultdict(collections.Counter)
changes = []
pairs = {}
for g in genomes:
    o, n = calls(OLD, g), calls(NEW, g)
    if o is None or n is None:
        tot["missing_report"] += 1
        continue
    f = meta.get(g, {}).get("FAMILY", "?")
    used = set()
    for a in o:
        m = next((i for i, b in enumerate(n) if i not in used and ov(a, b)), None)
        if m is None:
            k = "lost"
        else:
            used.add(m)
            b = n[m]
            k = "unchanged" if a["idiomorph"] == b["idiomorph"] else f"{a['idiomorph']}->{b['idiomorph']}"
            pairs[(g, a["contig"], a["start"])] = (a, b)
        tot[k] += 1
        fam[f][k] += 1
        if k != "unchanged":
            changes.append((g, f, meta.get(g, {}).get("SPECIES", ""), k, a["contig"], a["start"],
                            a.get("genes_found"), None if m is None else n[m].get("genes_found")))
    for i, b in enumerate(n):
        if i not in used:
            tot["new"] += 1
            fam[f]["new"] += 1
            changes.append((g, f, meta.get(g, {}).get("SPECIES", ""), "new", b["contig"], b["start"],
                            None, b.get("genes_found")))

print(f"genomes compared: {len(genomes)}")
print("overall:", dict(tot))
print("\nper family (only families with a change):")
for f, c in sorted(fam.items()):
    if any(k != "unchanged" for k in c):
        print(f"  {f:28s} {dict(c)}")

if SETS:
    rhizo = []
    for g in genomes:
        for a in calls(OLD, g) or []:
            gf = set(a["genes_found"])
            if a["idiomorph"] == "Plus" and "sexM" in gf and "sexP" not in gf:
                rhizo.append((g, a))
    c = collections.Counter()
    for g, a in rhizo:
        pr = pairs.get((g, a["contig"], a["start"]))
        c["lost" if pr is None else pr[1]["idiomorph"]] += 1
    print(f"\n33 Rhizopus sexM-only Plus calls ({len(rhizo)} found in OLD) -> NEW label: {dict(c)}")

    sexp_minus = {r["genome"] for r in csv.DictReader(open(CAND), delimiter="\t")
                  if r["status"] == "called" and r["idiomorph_label"] == "Minus" and r["tree_clade"] == "sexP"}
    print(f"\n{len(sexp_minus)} sexP-clade Minus genomes:")
    for g in sorted(sexp_minus):
        o, n = calls(OLD, g) or [], calls(NEW, g) or []
        sp = meta.get(g, {}).get("SPECIES", "")
        print(f"  {g[:30]:30s} {sp[:28]:28s} OLD {[x['idiomorph'] for x in o]} NEW {[x['idiomorph'] for x in n]}")
        for x in n:
            for e in x.get("idiomorph_resolutions") or []:
                if {e["winner"], e["loser"]} == {"sexM", "sexP"}:
                    print(f"      {e.get('basis')}: {e['winner']} {e.get('winner_model_score')} vs "
                          f"{e['loser']} {e.get('loser_model_score')} (id {e['winner_identity']:.1f}/{e['loser_identity']:.1f})")

    for fname in ("Lichtheimiaceae", "Syncephalastraceae"):
        gs = [g for g in genomes if meta.get(g, {}).get("FAMILY") == fname]
        co = sum(1 for g in gs if calls(OLD, g))
        cn = sum(1 for g in gs if calls(NEW, g))
        lab = collections.Counter(x["idiomorph"] for g in gs for x in calls(NEW, g) or [])
        print(f"\n{fname}: {len(gs)} genomes; called OLD {co}, NEW {cn}; NEW labels {dict(lab)}")

    rel = []
    for g in genomes:
        p = f"{NEW}/{g}/evidence_diagnostics.jsonl"
        if os.path.exists(p):
            rel += [json.loads(l) for l in open(p) if '"relaxed_call"' in l]
    reported = sum(1 for g in genomes for x in calls(NEW, g) or [] if x.get("detection_pass") == "relaxed")
    print(f"\nrelaxed-pass candidates {len(rel)}; reported {reported}; uncapped tier "
          f"{dict(collections.Counter(r['tier_uncapped'] for r in rel))}")
    for r in rel:
        print(f"  {r['genome_id'][:30]:30s} {r['contig']}:{r['start']} {r['idiomorph']} "
              f"genes={r['genes_found']} modelled={r['polished_genes']} uncapped={r['tier_uncapped']}")

print("\nall changes:")
for x in changes:
    print("  ", *x)
