"""The MAT-gene gate: is a classified core gene the MAT gene, or an HMG paralog?

Curator's ruling 2026-09-28 (option a), on the fragment leave-one-genus-out
validation (results/2026-09-28_validation_f3_f4/NOTE.md). The profile-HMM
idiomorph classifier TYPES a MAT gene reliably -- 0 wrong sexM/sexP calls in
540 held-out inputs -- but its margin does not tell a MAT gene from an HMG-box
paralog: at the 25-bit `min_margin`, 46 of 189 paralog proteins pass. Two
signals do separate them:

* the classifier's ABSOLUTE best score on a MODELLED protein: >= 100 bits for
  96/108 true MAT proteins but only 9/189 paralogs. On fragments (tblastn HSP
  translations) no score separates them;
* flank support: in the 293-genome Mucoromycota scan (f25cf70), every one of
  233 confident calls (model score >= 100, determined idiomorph) had at least
  one roster flank modelled at >= 45.8% identity, and 195 had two or more at
  >= 40%; the flanks of the paralog-like undetermined calls sat at 30.8-35.2%
  or came one at a time (a single glrA next to a Syncephalastrum paralog).

The rule, for a family whose roster names an `idiomorph_classifier`:

1. a model-typed call whose best classifier score is >= `mat_gene_min_score`
   (default 100 bits) is kept -- the gene is the MAT gene;
2. any other classified call (a lower model score, or a fragment- or
   mixed-typed call) is kept only with FLANK SUPPORT: at least
   `flank_support_min_genes` (default 2) DISTINCT roster flanking genes
   (role `flanking_*`) modelled at >= `flank_support_min_identity` (default
   40%) in the call;
3. otherwise the call is withheld with reason `mat_gene_gate`, and appears in
   `suppressed_loci` with the scores and the supporting flanks it had.

`min_margin` still decides TYPING; this gate decides only whether the call is
a MAT locus. Families without a classifier, calls with no classifier verdict,
and split-locus calls (which require roster flanks by construction) are left
alone. All three numbers are roster-overridable inside the
`idiomorph_classifier` block.
"""
from __future__ import annotations

from dataclasses import replace

from MATPredict.detect.polish import STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE

#: `DetectionResult.withheld_reason` for a call withheld by this gate.
WITHHELD_MAT_GENE_GATE = "mat_gene_gate"

DEFAULT_MAT_GENE_MIN_SCORE = 100.0
DEFAULT_FLANK_SUPPORT_MIN_IDENTITY = 40.0
DEFAULT_FLANK_SUPPORT_MIN_GENES = 2

_MODELLED_STATUSES = (STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE)


def _modelled(evidence) -> bool:
    """A polished model, or an annotated protein from the proteome fast path
    (as `flank_carried._modelled`)."""
    return evidence.status in _MODELLED_STATUSES or evidence.method == "diamond_proteome"


def supporting_flanks(result, min_identity: float) -> list[str]:
    """The distinct roster flanking genes modelled at >= `min_identity` in
    `result`, sorted."""
    return sorted({
        e.gene_name for e in result.gene_evidence
        if e.role.startswith("flanking") and _modelled(e)
        and e.identity is not None and e.identity >= min_identity
    })


def apply_mat_gene_gate(results, classifier_specs):
    """`(kept, withheld)` after the gate. `classifier_specs` maps each family
    key to its roster `idiomorph_classifier` block, or None."""
    kept, withheld = [], []
    for r in results:
        spec = classifier_specs.get(r.family_key)
        clf = r.idiomorph_classifier
        if not spec or not clf or r.split_locus:
            kept.append(r)
            continue
        min_score = float(spec.get("mat_gene_min_score", DEFAULT_MAT_GENE_MIN_SCORE))
        min_identity = float(spec.get("flank_support_min_identity",
                                      DEFAULT_FLANK_SUPPORT_MIN_IDENTITY))
        min_genes = int(spec.get("flank_support_min_genes", DEFAULT_FLANK_SUPPORT_MIN_GENES))
        best = max(clf.get("scores", {}).values(), default=0.0)
        if clf.get("classifier_input") == "model" and best >= min_score:
            kept.append(r)
            continue
        flanks = supporting_flanks(r, min_identity)
        if len(flanks) >= min_genes:
            kept.append(r)
            continue
        withheld.append(replace(
            r, withheld_reason=WITHHELD_MAT_GENE_GATE,
            withheld_detail={
                "best_score": round(best, 1),
                "classifier_input": clf.get("classifier_input"),
                "mat_gene_min_score": min_score,
                "supporting_flanks": flanks,
                "flank_support_min_genes": min_genes,
                "flank_support_min_identity": min_identity,
            },
        ))
    return kept, withheld
