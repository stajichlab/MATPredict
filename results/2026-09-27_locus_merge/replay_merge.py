"""Replay locus_merge.merge_overlapping on existing detection reports.

Usage: python replay_merge.py SRC_DIR DB_DIR REPORTS_GLOB...
Builds minimal DetectionResult objects from each report's `detected` list,
applies the real merge with the real roster groups, and prints merges.
"""
import glob, sys, collections
sys.path.insert(0, sys.argv[1])
import yaml
from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.locus_merge import merge_overlapping
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence

from pathlib import Path
fams = load_all_families(Path(sys.argv[2]))
groups = {f.key: (f.merge_group, f.merge_generic) for f in fams}
before = after = 0
merges = collections.Counter(); examples = []
for pat in sys.argv[3:]:
    for f in sorted(glob.glob(pat, recursive=True)):
        r = yaml.safe_load(open(f)) or {}
        res = []
        for x in r.get("detected") or []:
            ph, loc = x["family"].split(":", 1)
            ev = [GeneEvidence(gene_name=g["gene"], role=g.get("role", ""), contig=g["contig"],
                               start=g["start"], end=g["end"], strand=g.get("strand", "+"),
                               identity=g.get("identity") or 0, coverage=g.get("coverage"),
                               reference_record_id=g.get("reference_record", ""), method=g.get("method", ""))
                  for g in x.get("gene_evidence") or []]
            res.append(DetectionResult(family_key=FamilyKey(ph, loc), contig=x["contig"], start=x["start"],
                                       end=x["end"], confidence=x["confidence"], idiomorph=x["idiomorph"],
                                       ambiguous_with=[], genes_found=x.get("genes_found", []),
                                       genes_missing=x.get("genes_missing", []), fragmented=False,
                                       gene_evidence=ev, polished_genes=x.get("polished_genes", 0)))
        out = merge_overlapping(res, groups)
        before += len(res); after += len(out)
        for m in out:
            if m.merged_from:
                k = "+".join(sorted(d["family"].split(":")[1] for d in m.merged_from))
                merges[k] += 1
                if len(examples) < 6:
                    examples.append((f.split("/")[-2][:32], k, m.family_key.locus_name, m.confidence, m.idiomorph))
print("calls before", before, "after", after)
print("merges", dict(merges))
for e in examples: print("  ", e)
