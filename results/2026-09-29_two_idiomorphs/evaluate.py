"""Re-derive the two_idiomorphs statement from LCG reports (code bcd1e0d module).

Uses the committed module on the frozen-076afe4 LCG reports; calls are not
re-run, so no call can change. Genomes: LCG genomes/<Org>.sorted.fasta.
"""
import csv, glob, json, sys, yaml
from MATPredict.detect.assembly_gap import _contig_sequences
from MATPredict.detect.two_idiomorphs import calls_from_report, two_idiomorph_statements

LCG = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-28_lcg_holdout"
GEN = "/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/genomes"
FAM = {"Mucoromycota:MAT"}
LIT = {"Syzygites_megalocarpus_SC16", "Syzygites_sp._MES_3091",
       "Zygorhynchus_heterogamus_Vuillemin_NRRL_1489",
       "Zygorhynchus_moelleri_NRRL_1498", "Zygorhynchus_moelleri_NRRL_1625",
       "Zygorhynchus_moelleri_NRRL_3138", "Zygorhynchus_moelleri_Vuillemin_NRRL_1497"}
# conventionally described as heterothallic (check group; general literature)
HET = ("Mucor_hiemalis", "Mucor_racemosus", "Rhizopus_nigricans", "Rhizopus_oligosporus",
       "Rhizopus_chinensis", "Rhizopus_rhizopodiformis", "Backusella", "Circinella",
       "Mucor_indicus", "Mucor_rouxii", "Mucor_rouxianus")
HOMO_REPORTED = ("Mucor_genevensis", "Rhizopus_azygosporus", "Mucor_azygosporus",
                 "Zygorhynchus", "Syzygites")

def group(org):
    if any(org.startswith(p) for p in HOMO_REPORTED): return "reported_homothallic"
    if any(org.startswith(p) for p in HET): return "conventionally_heterothallic"
    return "unknown"

rows, stmts = [], {}
orgs = sorted({f.split("/")[-2] for f in glob.glob(f"{LCG}/runs2/*/detection_report.yaml")})
for org in orgs:
    if org.startswith("Mycotypha_africana_NRRL_2978"):
        continue  # training strain, excluded
    doc = yaml.safe_load(open(f"{LCG}/runs2/{org}/detection_report.yaml"))
    calls = calls_from_report(doc)
    ids = {c["idiomorph"] for c in calls if c["family"] in FAM} - {None, "undetermined"}
    if not ({"Plus", "Minus"} <= ids or org in LIT):
        continue
    wanted = {c["contig"] for c in calls} | {g["contig"] for c in calls for g in c["gene_evidence"]}
    seqs = _contig_sequences(f"{GEN}/{org}.sorted.fasta", wanted)
    st = two_idiomorph_statements(calls, FAM, seqs, genetic_code=doc.get("genetic_code") or 1)
    stmts[org] = {"statements": st, "calls": [(c["idiomorph"], c["confidence"], c["contig"], c["start"], c["end"], c["margin"]) for c in calls],
                  "suppressed_near": [ {k: s.get(k) for k in ("contig","start","end","idiomorph","withheld_reason")}
                                       for s in doc.get("suppressed_loci") or [] if s.get("family") in FAM and s.get("idiomorph") not in (None,"undetermined")]}
    s = st[0] if st else None
    sf = (s or {}).get("evidence", {}).get("shared_flanks", []) if s else []
    rows.append({
        "org": org, "group": group(org), "literature_positive": org in LIT,
        "idiomorphs_called": "+".join(sorted(ids)),
        "arrangement": s["arrangement"] if s else "",
        "plus_conf": s["plus_call"]["confidence"] if s else "", "minus_conf": s["minus_call"]["confidence"] if s else "",
        "plus_margin": s["plus_call"]["margin"] if s else "", "minus_margin": s["minus_call"]["margin"] if s else "",
        "plus_input": s["plus_call"]["classifier_input"] if s else "", "minus_input": s["minus_call"]["classifier_input"] if s else "",
        "plus_flanks": ",".join(f["gene"] for f in s["plus_call"]["flanks"]) if s else "",
        "minus_flanks": ",".join(f["gene"] for f in s["minus_call"]["flanks"]) if s else "",
        "shared_flank_identity": ";".join(f"{x['gene']}={x['protein_identity']}" for x in sf),
        "gc_diff": s["evidence"].get("gc_difference_pct") if s else "",
        "calls_per_idiomorph": json.dumps(s["evidence"]["calls_per_idiomorph"]) if s else "",
        "weak_calls": ",".join(s["evidence"]["weak_calls"]) if s else "",
        "supported_causes": ",".join(s["supported_causes"]) if s else "",
    })
with open("per_genome.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
json.dump(stmts, open("statements.json", "w"), indent=1, default=str)
print(len(rows), "rows")
