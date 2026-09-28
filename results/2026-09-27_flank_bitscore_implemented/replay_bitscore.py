"""Replay the 39-bit flank-carried rule through the REAL function (code and db
of the frozen worktree run-8d80bed) on the audit's 108 changed calls
(changed.tsv). Bitscores are the audit's own tblastn values (classified3.tsv
here_bits, best in-locus per gene) -- the cap6 reports predate
GeneEvidence.bitscore. Compares with the E <= 1e-5 rule on the same inputs.
usage (from the audit dir): PYTHONPATH=<run-8d80bed>/src python replay_bitscore.py"""
import collections, csv, os
from pathlib import Path

from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.flank_carried import apply_flank_carried_rule
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence

DB = Path(os.environ["RUN_DB"])
fams = load_all_families(DB)
windows = {f.key: f.flank_carried_window_bp for f in fams}
floors = {f.key: f.flank_carried_min_bitscore for f in fams}

bits, ev = collections.defaultdict(dict), collections.defaultdict(dict)
for r in csv.DictReader(open("classified3.tsv"), delimiter="\t"):
    k = (r["panel"], r["genome"], r["mat_family"], r["call"])
    if r["here_bits"]:
        bits[k][r["gene"]] = float(r["here_bits"])
    if r["here_evalue"]:
        ev[k][r["gene"]] = float(r["here_evalue"])

def evidence(detail, role, contig, b):
    out = []
    for d in filter(None, detail.split(";")):
        gene, pos, ident, cov, status = d.split(":")
        s, e = map(int, pos.split("-"))
        out.append(GeneEvidence(gene, role, contig, s, e, "+", float(ident),
                                None if cov == "None" else float(cov), "rec", "tblastn_genome",
                                status=status, bitscore=b.get(gene) if role == "core_MAT" else None))
    return out

def e_rule(res, k, window):
    """The E <= 1e-5 rule on the same evidence (strongest = lowest E)."""
    core = [e for e in res.gene_evidence if e.role == "core_MAT"]
    fl = [e for e in res.gene_evidence if e.role.startswith("flanking")]
    lo = min(e.start for e in fl); hi = max(e.end for e in fl)
    if not core: return "withheld"
    s = min(core, key=lambda e: ev[k].get(e.gene_name, 9e9))
    d = max(lo - s.end, s.start - hi, 0)
    return "low" if ev[k].get(s.gene_name, 9e9) <= 1e-5 and d <= window else "withheld"

tally = collections.Counter(); rows = []
for r in csv.DictReader(open("changed.tsv"), delimiter="\t"):
    call = f'{r["contig"]}:{r["start"]}-{r["end"]}'
    phylum, locus = r["mat_family"].split(":")
    key = FamilyKey(phylum, locus); k = (r["panel"], r["genome"], r["mat_family"], call)
    gene_ev = (evidence(r["core_detail"], "core_MAT", r["contig"], bits[k])
               + evidence(r["flank_detail"], "flanking_conserved", r["contig"], {}))
    res = DetectionResult(key, r["contig"], int(r["start"]), int(r["end"]), r["confidence"],
                          r["idiomorph"], [], [], [], False, gene_evidence=gene_ev)
    kept, _ = apply_flank_carried_rule([res], windows, floors)
    new = "low" if kept else "withheld"
    old = e_rule(res, k, windows.get(key, 3000))
    grp = "Ascomycota" if r["panel"] == "cap6" else "Mucoromycota-group"
    tally[(grp, "bits", new)] += 1; tally[(grp, "E", old)] += 1
    rows.append((grp, r["genome"], r["species"], r["mat_family"], call, old, new))

for grp in ("Ascomycota", "Mucoromycota-group"):
    print(f"{grp}: 39-bit rule kept {tally[(grp,'bits','low')]}, withheld {tally[(grp,'bits','withheld')]}"
          f" | E<=1e-5 rule kept {tally[(grp,'E','low')]}, withheld {tally[(grp,'E','withheld')]}")
print("\nnewly kept by the 39-bit rule (withheld under E):")
for x in rows:
    if x[5] == "withheld" and x[6] == "low": print("  ", x)
print("\nnewly withheld by the 39-bit rule (kept under E):")
for x in rows:
    if x[5] == "low" and x[6] == "withheld": print("  ", x)
