"""Compare Wallemiales calls: one-record run (a6770ca) vs two-record run (3d8a755).

Version per genome (v1 = CBS 633.66 type, other) comes from the independent
check in results/2026-09-27_puccinio_followup/wallemia/versions.tsv.
"""
import collections, csv, glob, os, yaml

BASE = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
OLD = f"{BASE}/2026-09-27_puccinio_followup/stepB_wallemia/runs"
NEW = f"{BASE}/2026-09-27_wallemia_allele2/run_3d8a755/runs"
ver = {r["asmid"]: r for r in csv.DictReader(open(f"{BASE}/2026-09-27_puccinio_followup/wallemia/versions.tsv"), delimiter="\t")}


def calls(d, a):
    f = f"{d}/{a}/detection_report.yaml"
    if not os.path.exists(f):
        return None
    r = yaml.safe_load(open(f))
    out = []
    for x in r.get("detected") or []:
        if x["family"].endswith("wallMAT"):
            ge = {g["gene"]: (round(g["identity"], 1) if g.get("identity") is not None else None, g["status"])
                  for g in x.get("gene_evidence", [])}
            out.append(dict(idio=x["idiomorph"], conf=x["confidence"], cls=x["locus_class"],
                            cands=x.get("idiomorph_candidates"), genes=ge, contig=x["contig"], start=x["start"]))
    return out


rows = []
for a, v in sorted(ver.items(), key=lambda kv: (kv[1]["species"], kv[1]["version"])):
    o, n = calls(OLD, a), calls(NEW, a)
    vv = "v1" if v["version"].startswith("v1") else "other"
    rows.append((a, v["species"], v["strain"], vv, o, n))

tally = collections.Counter()
with open(f"{BASE}/2026-09-27_wallemia_allele2/per_genome.tsv", "w") as fo:
    fo.write("asmid\tspecies\tstrain\tversion_check\told_call\tnew_call\tnew_idiomorph_scores\tnew_STE3\tnew_STE3v2\n")
    for a, sp, st, vv, o, n in rows:
        oc = ";".join(f"{c['idio']}/{c['conf']}/{c['cls']}" for c in (o or [])) or "none"
        nc = ";".join(f"{c['idio']}/{c['conf']}/{c['cls']}" for c in (n or [])) or ("NO_REPORT" if n is None else "none")
        sc = ";".join(str(c["cands"]) for c in (n or []))
        s1 = ";".join(str(c["genes"].get("STE3")) for c in (n or []))
        s2 = ";".join(str(c["genes"].get("STE3v2")) for c in (n or []))
        fo.write(f"{a}\t{sp}\t{st}\t{vv}\t{oc}\t{nc}\t{sc}\t{s1}\t{s2}\n")
        idios = sorted({c["idio"] for c in (n or [])})
        tally[(vv, "+".join(idios) or "none", "+".join(sorted({c['conf'] for c in (n or [])})) or "-")] += 1
        tally[(sp, vv, "+".join(idios) or "none")] += 0
    for k, c in sorted(tally.items(), key=str):
        if c:
            print(k, c)
sp_split = collections.Counter()
for a, sp, st, vv, o, n in rows:
    idios = "+".join(sorted({c["idio"] for c in (n or [])})) or "none"
    sp_split[(sp, idios)] += 1
print()
for k, c in sorted(sp_split.items()):
    print(k, c)
changed_v1 = [(a, sp) for a, sp, st, vv, o, n in rows if vv == "v1" and
              sorted((c['conf'], c['cls']) for c in (o or [])) != sorted((c['conf'], c['cls']) for c in (n or []))]
print("\nv1-type genomes whose confidence/class changed:", len(changed_v1), changed_v1[:10])
