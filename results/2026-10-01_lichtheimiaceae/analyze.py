"""Steps 1, 3, 4: typed sexM/sexP candidates per genome, Plus/Minus by species, H1/H2 replay.

Reads all.tsv + out/<gid>/{loci,clf_models,clf_annot,mp}.tsv and detection reports on current
code (LCG: 2026-10-01_heldout_rerun/lcg_runs at 52b3ff9; BFD: 2026-09-30_cap_v3 cand runs at
cf5cb1e = V3; Jena: 2026-10-01_heldout_rerun/jena_scaffolds at 52b3ff9).
Writes candidates.tsv, per_genome.tsv, by_genus.tsv, h1_h2_gains.tsv, summary.txt.
"""
import csv, yaml
from pathlib import Path
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent
R = HERE.parent
GATE, MARGIN = 98.4, 25.0
REP = {"LCG": R / "2026-10-01_heldout_rerun/lcg_runs", "BFD": R / "2026-09-30_cap_v3/regression/cand/mucoromycota/runs",
       "Jena": R / "2026-10-01_heldout_rerun/jena_scaffolds"}


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t")) if Path(p).exists() else []


def typed(m, p, p1):
    top = max(m, p)
    if p1 >= top + MARGIN:
        return "paralog_P1"
    if abs(m - p) < MARGIN:
        return "undetermined"
    lab = "Minus" if m > p else "Plus"
    return lab if top >= GATE else lab + "_belowgate"


def report(setn, gid):
    rid = gid[5:] if setn == "Jena" else gid
    f = REP[setn] / rid / "detection_report.yaml"
    return yaml.safe_load(open(f)) if f.exists() else None


rows = tsv(HERE / "all.tsv")
famg = {x["genome_id"]: x["genus"] for x in tsv(HERE / "genomes.tsv")}
for r in rows:
    if r["genome_id"] in famg:
        r["genus"] = famg[r["genome_id"]]
        r["note"] = "lichtheimiaceae"
cand_out, gen_out, gains = [], [], []
for r in rows:
    gid, setn = r["genome_id"], r["set"]
    od = HERE / "out" / gid
    if not (od / "done").exists():
        continue
    loci = {x["mid"]: x for x in tsv(od / "loci.tsv")}
    best_locus = max(loci.values(), key=lambda x: float(x["top_tblastn_bits"] or 0), default=None)
    cands = []
    for c in tsv(od / "clf_models.tsv"):
        m, p, p1 = float(c["sexM"]), float(c["sexP"]), float(c["P1"])
        L = loci.get(c["mid"], {})
        cands.append(dict(genome_id=gid, set=setn, source="model", id=c["mid"], contig=L.get("contig"), start=L.get("start"),
                          end=L.get("end"), sexM=m, sexP=p, P1=p1, call=typed(m, p, p1), model_len=L.get("model_len"),
                          tblastn_bits=L.get("top_tblastn_bits"), best_locus=bool(best_locus and L is best_locus)))
    for c in tsv(od / "clf_annot.tsv"):
        m, p, p1 = float(c["sexM"]), float(c["sexP"]), float(c["P1"])
        if max(m, p) < 50:
            continue
        cands.append(dict(genome_id=gid, set=setn, source="annot", id=c["protein"], contig="", start="", end="", sexM=m, sexP=p,
                          P1=p1, call=typed(m, p, p1), model_len="", tblastn_bits="", best_locus=False))
    cand_out += cands
    strong = sorted({c["call"] for c in cands if c["call"] in ("Plus", "Minus")})
    strong_models = [c for c in cands if c["source"] == "model" and c["call"] in ("Plus", "Minus")]
    rep = report(setn, gid)
    calls = [(d["contig"], int(d["start"]), int(d["end"]), d.get("idiomorph")) for d in (rep or {}).get("detected") or []]
    called = sorted({c[3] for c in calls}) if calls else []
    # H2: best tblastn locus is core-only/uncalled, strong typed model there, no call overlapping it
    h2 = ""
    if best_locus:
        bm = [c for c in strong_models if c["best_locus"]]
        overl = any(c[0] == best_locus["contig"] and c[1] <= int(best_locus["end"]) + 5000 and c[2] >= int(best_locus["start"]) - 5000 for c in calls)
        if bm and not overl:
            h2 = bm[0]["call"]
            gains.append(dict(hyp="H2", genome_id=gid, set=setn, genus=r["genus"], species=r["species"], lichtheimiaceae=r["note"],
                              call=h2, sexM=bm[0]["sexM"], sexP=bm[0]["sexP"], P1=bm[0]["P1"], contig=best_locus["contig"],
                              start=best_locus["start"], already_called=";".join(called), detail="best tblastn locus not called"))
    # H1: gate-withheld loci with a strong rnhA miniprot hit within 20 kb
    h1 = []
    rn = [x for x in tsv(od / "mp.tsv") if x["gene"] == "rnhA" and float(x["identity"] or 0) >= 0.5]
    for s in (rep or {}).get("suppressed_loci") or []:
        if s.get("withheld_reason") != "mat_gene_gate":
            continue
        clf = s.get("idiomorph_classifier") or {}
        if clf.get("idiomorph") not in ("Plus", "Minus"):
            continue
        near = [x for x in rn if x["contig"] == s["contig"] and int(x["start"]) <= s["end"] + 20000 and int(x["end"]) >= s["start"] - 20000]
        if near:
            idn = max(float(x["identity"]) for x in near)
            h1.append(clf["idiomorph"])
            gains.append(dict(hyp="H1", genome_id=gid, set=setn, genus=r["genus"], species=r["species"], lichtheimiaceae=r["note"],
                              call=clf["idiomorph"], sexM=clf.get("scores", {}).get("Minus"), sexP=clf.get("scores", {}).get("Plus"),
                              P1=(clf.get("paralog_scores") or {}).get("P1"), contig=s["contig"], start=s["start"],
                              already_called=";".join(called), detail=f"rnhA id {idn:.2f}; input {clf.get('classifier_input')}"))
    gen_out.append(dict(genome_id=gid, set=setn, genus=r["genus"], species=r["species"], lichtheimiaceae=r["note"],
                        label=r["label"], label_disputed=r["label_disputed"], typed_strong=";".join(strong) or "none",
                        reported_calls=";".join(called) or "none", H2=h2, H1=";".join(h1)))

for name, data in (("candidates.tsv", cand_out), ("per_genome.tsv", gen_out), ("h1_h2_gains.tsv", gains)):
    if data:
        with open(HERE / name, "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=list(data[0]), delimiter="\t"); w.writeheader(); w.writerows(data)

# by genus (Lichtheimiaceae only)
lic = [g for g in gen_out if g["lichtheimiaceae"]]
bg = defaultdict(Counter)
for g in lic:
    k = g["genus"]
    bg[k]["genomes"] += 1
    bg[k]["typed_" + g["typed_strong"]] += 1
    bg[k]["reported_" + g["reported_calls"]] += 1
    if g["H2"]: bg[k]["H2_gain"] += 1
    if g["H1"]: bg[k]["H1_gain"] += 1
with open(HERE / "by_genus.tsv", "w") as fo:
    keys = sorted({k for c in bg.values() for k in c})
    fo.write("genus\t" + "\t".join(keys) + "\n")
    for gname, c in sorted(bg.items()):
        fo.write(gname + "\t" + "\t".join(str(c.get(k, 0)) for k in keys) + "\n")
with open(HERE / "summary.txt", "w") as fo:
    fo.write(f"genomes analysed: {len(gen_out)} (Lichtheimiaceae {len(lic)})\n")
    for hyp in ("H1", "H2"):
        gl = [x for x in gains if x["hyp"] == hyp]
        newg = {x["genome_id"] for x in gl if not x["already_called"]}
        fo.write(f"{hyp}: loci {len(gl)}; genomes with no current call gained {len(newg)}; "
                 f"in Lichtheimiaceae {sum(1 for x in gl if x['lichtheimiaceae'])}; outside {sum(1 for x in gl if not x['lichtheimiaceae'])}\n")
        fo.write(f"   by set {Counter(x['set'] for x in gl)}; calls {Counter(x['call'] for x in gl)}\n")
print(open(HERE / "summary.txt").read())
