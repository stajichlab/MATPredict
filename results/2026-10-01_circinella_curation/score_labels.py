"""Score Circinella-group LCG strains against the curator's label treatment
(results/2026-10-01_circinella_label_tree/label_treatment.tsv), base vs cand.

Truth = the strain's tree clade (sexM -> Minus, sexP -> Plus), which is what
label_treatment.tsv rules for every strain (label reliable, wrong or
inverted). PROVISIONAL strains are listed but NOT scored. Also lists every
Circinella-group LCG genome (by genus) with its base and cand MAT calls.

Usage: python score_labels.py LABELS LCG_DIR STRAIN_TABLE > label_score.tsv
"""
import csv, sys
from pathlib import Path
import yaml

labels, lcg, strains = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
GENERA = ("Circinella", "Thamnostylum", "Fennellomyces", "Zychaea")
CLADE = {"sexM": "Minus", "sexP": "Plus"}


def calls(side, org):
    p = lcg / side / org / "detection_report.yaml"
    if not p.exists():
        return "no_report"
    rep = yaml.safe_load(p.read_text()) or {}
    out = [f"{d.get('idiomorph')}:{d.get('confidence')}" + ("*grp" if d.get("mat_gene_gate_group") else "")
           for d in rep.get("detected") or [] if str(d.get("family", "")).endswith(":MAT")]
    held = [f"withheld:{d.get('withheld_reason')}" for d in rep.get("suppressed_loci") or []
            if str(d.get("family", "")).endswith(":MAT") and d.get("withheld_reason") == "mat_gene_gate"]
    return ";".join(out) if out else ("none" + ("(" + ";".join(sorted(set(held))) + ")" if held else ""))


def verdict(call, truth):
    ids = {c.split(":")[0] for c in call.split(";") if ":" in c and not c.startswith("withheld")}
    if not ids or call.startswith("none") or call == "no_report":
        return "uncalled"
    if ids == {"undetermined"}:
        return "undetermined"
    if ids == {truth}:
        return "correct"
    if truth in ids:
        return "both"
    return "wrong"


print("strain\ttruth\tprovisional\tbase\tbase_verdict\tcand\tcand_verdict")
tot = {"base": {}, "cand": {}}
for r in csv.DictReader(open(labels), delimiter="\t"):
    truth = CLADE[r["tree_clade"]]
    prov = "PROVISIONAL" in r["treatment"]
    row = [r["strain"], truth, "yes" if prov else "no"]
    for side in ("base", "cand"):
        c = calls(side, r["strain"])
        v = verdict(c, truth)
        row += [c, v]
        if not prov:
            tot[side][v] = tot[side].get(v, 0) + 1
    print(*row, sep="\t")
print("#scored (non-provisional)\tbase=" + str(tot["base"]) + "\tcand=" + str(tot["cand"]))
print("#all Circinella-group LCG genomes (not scored)")
print("genome\tstrain_table_type\tbase\tcand")
types = {r["strain"]: r["rnhA_gene"] for r in csv.DictReader(open(strains), delimiter="\t")}
for d in sorted((lcg / "cand").iterdir()):
    if d.name.split("_")[0] in GENERA:
        print(d.name, types.get(d.name, ""), calls("base", d.name), calls("cand", d.name), sep="\t")
