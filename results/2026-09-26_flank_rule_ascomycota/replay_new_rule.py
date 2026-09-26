"""Replay the REVISED flank-carried rule (src/MATPredict/detect/flank_carried.py)
through the real function, with the curated per-family windows from db/, on the
108 calls the first rule changed (changed.tsv). E-values come from this audit's
tblastn (classified3.tsv, best in-locus e-value per gene), because the reports
predate GeneEvidence.evalue. Run from this directory with PYTHONPATH=<repo>/src.
usage: python replay_new_rule.py > replay_new_rule.txt"""
import collections, csv
from pathlib import Path

from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.flank_carried import apply_flank_carried_rule
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence

DB = Path(__file__).resolve().parents[2] / "db"
windows = {f.key: f.flank_carried_window_bp for f in load_all_families(DB)}

ev = collections.defaultdict(dict)
for r in csv.DictReader(open("classified3.tsv"), delimiter="\t"):
    if r["here_evalue"]:
        ev[(r["panel"], r["genome"], r["mat_family"], r["call"])][r["gene"]] = float(r["here_evalue"])

def evidence(detail, role, contig, evalues):
    out = []
    for d in filter(None, detail.split(";")):
        gene, pos, ident, cov, status = d.split(":")
        s, e = map(int, pos.split("-"))
        out.append(GeneEvidence(gene, role, contig, s, e, "+", float(ident),
                                None if cov == "None" else float(cov), "rec", "tblastn_genome",
                                status=status, evalue=evalues.get(gene) if role == "core_MAT" else None))
    return out

tally = collections.Counter()
rows = []
for r in csv.DictReader(open("changed.tsv"), delimiter="\t"):
    call = f'{r["contig"]}:{r["start"]}-{r["end"]}'
    phylum, locus = r["mat_family"].split(":")
    key = FamilyKey(phylum, locus)
    evalues = ev[(r["panel"], r["genome"], r["mat_family"], call)]
    gene_ev = (evidence(r["core_detail"], "core_MAT", r["contig"], evalues)
               + evidence(r["flank_detail"], "flanking_conserved", r["contig"], {}))
    res = DetectionResult(key, r["contig"], int(r["start"]), int(r["end"]), r["confidence"],
                          r["idiomorph"], [], [], [], False, gene_evidence=gene_ev)
    kept, withheld = apply_flank_carried_rule([res], windows)
    outcome = "low" if kept else "withheld"
    grp = "Ascomycota" if r["panel"] == "cap6" else "Mucoromycota-group"
    tally[(grp, outcome)] += 1
    rows.append((grp, r["genome"], r["species"], r["mat_family"], call, windows.get(key), outcome))

sim = {}
for s in csv.DictReader(open("simulated.tsv"), delimiter="\t"):
    sim[(s["genome"], s["mat_family"], s["call"])] = s["B_20kb_E1e-5"]
diff = [x for x in rows if sim.get((x[1], x[3], x[4])) != x[6]]
for grp in ("Ascomycota", "Mucoromycota-group"):
    print(f"{grp}: low {tally[(grp, 'low')]}, withheld {tally[(grp, 'withheld')]}")
print("kept real loci:", [(x[2], x[3], x[6]) for x in rows
                          if x[2] in ("Mucor irregularis", "Trigonopsis variabilis")])
print(f"differences from simulated.tsv B_20kb_E1e-5 (window 20 kb everywhere): {len(diff)}")
for x in diff:
    print("  ", x, "simulated:", sim.get((x[1], x[3], x[4])))
