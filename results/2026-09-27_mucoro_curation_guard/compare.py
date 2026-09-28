"""Compare the guard+Syncephalastrum run (2e9aa97) against the split-locus run (4174440)."""
import csv, glob, collections, yaml, sys
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
OLD = f"{R}/2026-09-27_split_locus/Mucoromycota_4174440/runs"
NEW = f"{R}/2026-09-27_mucoro_curation_guard/Mucoromycota_2e9aa97/runs"
samp = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
def load(d):
    out = {}
    for f in glob.glob(f"{d}/*/detection_report.yaml"):
        out[f.split("/")[-2]] = yaml.safe_load(open(f))
    return out
old, new = load(OLD), load(NEW)
def calls(r):
    return sorted((x["contig"], x["start"], x["idiomorph"], x["confidence"], x["locus_class"],
                   (x.get("idiomorph_classifier") or {}).get("margin")) for x in (r.get("detected") or []))
fam_of = lambda a: samp.get(a, {}).get("FAMILY") or "?"
tally = collections.Counter(); lines = []
guard = collections.Counter(); guard_list = []
for a in sorted(set(old) | set(new)):
    o, n = calls(old.get(a, {})), calls(new.get(a, {}))
    f = fam_of(a)
    grp = f if f in ("Umbelopsidaceae", "Lichtheimiaceae", "Syncephalastraceae") else "other"
    for x in (new.get(a, {}).get("suppressed_loci") or []):
        if x.get("withheld_reason") == "secondary_call_classifier_undetermined":
            guard[grp] += 1; guard_list.append((a, samp.get(a, {}).get("SPECIES"), x["contig"], x["start"]))
    ok = {(c[0], c[1]) for c in o}; nk = {(c[0], c[1]) for c in n}
    oi = collections.Counter(c[2] for c in o); ni = collections.Counter(c[2] for c in n)
    if not o and n: tally[(grp, "gained_genome")] += 1
    if o and not n: tally[(grp, "lost_genome")] += 1
    tally[(grp, "calls_old")] += len(o); tally[(grp, "calls_new")] += len(n)
    if o != n:
        lines.append(f"{grp:18s} {samp.get(a,{}).get('SPECIES')} {samp.get(a,{}).get('STRAIN')} {a}\n   old: {o}\n   new: {n}")
print("per group:")
for g in ("Umbelopsidaceae", "Lichtheimiaceae", "Syncephalastraceae", "other"):
    print(f"  {g:18s} calls {tally[(g,'calls_old')]} -> {tally[(g,'calls_new')]}; genomes gained {tally[(g,'gained_genome')]}, lost {tally[(g,'lost_genome')]}; guard-withheld {guard[g]}")
print("genomes with any call:", sum(1 for a in old if old[a].get('detected')), "->", sum(1 for a in new if new[a].get('detected')))
print("\nguard-withheld calls:")
for x in guard_list: print("  ", *x)
print("\nchanged genomes:")
print("\n".join(lines))
