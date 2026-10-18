"""Both-genes homothallic rule: a COMPLETE alpha-box gene and a COMPLETE HMG-box gene at one locus.

Behind `--homothallic-policy both_genes` (default `off`: nothing here runs and the reports are unchanged).
Curator approval 2026-10-18 (question 16). Background and measurements:
analysis/2026-10-18_dothideomycetes-mat-audit.md and analysis/2026-10-18_homothallic-loci-patch.md.

Why a second rule next to `idiomorph.classify_locus`. The existing different-domains path needs both genes to be
`full_length_models`: modelled over >= 50% of their REFERENCE protein, and in no idiomorph cross-match resolution.
On the 298 Dothideomycetes genomes whose locus holds both genes but was reported as one idiomorph it failed for
four independent reasons: the relaxed pass calls `classify_locus` without `full_length_models` at all (12); the HMG
gene lost an overlap resolution to the optional MAT1-1-3 or was never listed (97); the MAT1-2-1/MAT1-1-3 cross-match
exclusion (188, the largest); and the model covers under 50% of a distant reference (169, overlapping the others).

This rule does not use the reference at all. A gene counts when
  * its model is polished (agree, single or disagree), whatever resolution it took part in or lost;
  * its translated ORF is complete: the longest stop-free stretch is >= MIN_ORF_FRACTION of the model's coding
    length and >= MIN_PROTEIN_AA residues;
  * the shipped gene-identity profiles confirm WHICH gene it is: its own profile scores >= MIN_PROFILE_BITS and beats
    the other profile by >= PROFILE_MIN_MARGIN. An HMG protein that is not the MAT1-2-1 profile's (a paralog, or the
    HMG gene of MAT1-1-3) does not pass; the HMG floor sits above every curated MAT1-1-3 and below every curated
    MAT1-2-1 (leave-one-genus-out), see CALIBRATION below.
and the two genes are different models on one contig, not overlapping, within the family's
`max_homothallic_separation_bp`. The class stays a CANDIDATE: homothallism is never asserted and `idiomorph` keeps
the vote's value (see the note).

CALIBRATION (curated records only, no BFD genome; results/2026-10-18_homothallic_patch/calibration.txt):
  MAT1-1-1 on the alpha profile, leave-one-genus-out, 22 records: 102.4 to 263.4 bits (cross score <= 18).
  MAT1-2-1 on the HMG profile, leave-one-genus-out, 26 records: 120.5 to 291.5 bits (cross score <= 15).
  Curated MAT1-1-3 (HMG-containing, MAT1-1 idiomorph; 4 proteins, 162-324 aa) on the HMG profile: 98.0, 100.9, 100.9, 51.2.
  So MIN_PROFILE_BITS = {alpha: 100, HMG: 110}.  Length floor 200 aa: 23 of 26 curated MAT1-1-1 are >= 296 aa and 24 of 29 MAT1-2-1 are
  >= 223 aa; the shorter ones are partial deposits.
All of these are PROVISIONAL, like the other curator-tunable thresholds.
"""
from __future__ import annotations

from dataclasses import dataclass, field

POLICY_OFF = "off"
POLICY_BOTH_GENES = "both_genes"
HOMOTHALLIC_POLICIES = (POLICY_OFF, POLICY_BOTH_GENES)
DEFAULT_HOMOTHALLIC_POLICY = POLICY_OFF

#: Profile floors by gene class (bits), from the curated records only.
MIN_PROFILE_BITS = {"alpha_box": 100.0, "HMG_box": 110.0}
#: A gene's own profile must beat the other profile by this many bits (the same 25 as the partial-core rule).
PROFILE_MIN_MARGIN = 25.0
#: Longest stop-free stretch of the translated model over its coding length.
MIN_ORF_FRACTION = 0.9
#: Complete ORF at least this long (aa).
MIN_PROTEIN_AA = 200
MODELLED_STATUSES = ("polished_agree", "polished_single", "polished_disagree")


@dataclass(frozen=True)
class CoreModel:
    """One modelled core gene of a locus, with what the rule needs to judge it."""

    gene: str
    contig: str
    start: int
    end: int
    protein_aa: int
    """Longest stop-free stretch of the translated model."""
    cds_aa: int
    """Coding length of the model's exons / 3."""
    scores: dict[str, float] = field(default_factory=dict)
    """Best bits on each gene-class profile, keyed `alpha_box` / `HMG_box`."""
    lost_to: str | None = None


def confirmed_class(m: CoreModel) -> str | None:
    """`alpha_box` or `HMG_box` when this model is a complete, profile-confirmed gene of that class, else None."""
    if m.cds_aa <= 0 or m.protein_aa < MIN_PROTEIN_AA or m.protein_aa < MIN_ORF_FRACTION * m.cds_aa:
        return None
    for cls, other in (("alpha_box", "HMG_box"), ("HMG_box", "alpha_box")):
        own, cross = m.scores.get(cls, 0.0), m.scores.get(other, 0.0)
        if own >= MIN_PROFILE_BITS[cls] and own - cross >= PROFILE_MIN_MARGIN:
            return cls
    return None


def evaluate(models: list[CoreModel], max_separation_bp: int) -> dict | None:
    """The note of a both-genes candidate, or None. Picks the pair with the highest summed profile score."""
    alpha = [m for m in models if confirmed_class(m) == "alpha_box"]
    hmg = [m for m in models if confirmed_class(m) == "HMG_box"]
    best = None
    for a in alpha:
        for h in hmg:
            if a.contig != h.contig:
                continue
            separation = max(a.start, h.start) - min(a.end, h.end)
            if separation <= 0 or separation > max_separation_bp:
                continue  # overlapping models are one region seen twice; far apart is not one locus
            key = a.scores["alpha_box"] + h.scores["HMG_box"]
            if best is None or key > best[0]:
                best = (key, a, h, separation)
    if best is None:
        return None
    _, a, h, separation = best
    return {
        "decision": "homothallic_candidate",
        "policy": POLICY_BOTH_GENES,
        "alpha_box_gene": {"gene": a.gene, "protein_aa": a.protein_aa, "alpha_bits": round(a.scores["alpha_box"], 1),
                           "hmg_bits": round(a.scores.get("HMG_box", 0.0), 1), "lost_to": a.lost_to},
        "hmg_box_gene": {"gene": h.gene, "protein_aa": h.protein_aa, "hmg_bits": round(h.scores["HMG_box"], 1),
                         "alpha_bits": round(h.scores.get("alpha_box", 0.0), 1), "lost_to": h.lost_to},
        "separation_bp": separation,
        "reason": (f"complete alpha-box protein ({a.gene}, {a.protein_aa} aa, {a.scores['alpha_box']:.0f} bits) and complete "
                   f"HMG-box protein ({h.gene}, {h.protein_aa} aa, {h.scores['HMG_box']:.0f} bits) on {a.contig}, "
                   f"{separation} bp apart; a candidate, not a statement of homothallism"),
    }


def note_for_evidence(evidence, family, clf, translate, max_separation_bp: int) -> dict | None:
    """`evaluate` over a locus's gene evidence.

    `evidence` are the call's `GeneEvidence`; `clf` the gene-identity profiles (`partial_core.load_core_profiles`);
    `translate(span)` the protein a model encodes (`span` carries contig, start, end, strand, exons), so the caller's
    genome access and genetic code stay outside. Every polished core gene is judged, including one that lost an
    idiomorph cross-match resolution: whether it is the MAT1-2-1 HMG gene is decided by the profile, not by the name.
    """
    import types

    from MATPredict.detect.classifier import score_proteins

    class_of = {idiomorph: next((g.get("gene_class") for g in family.genes if g["name"] == gene), None)
                for idiomorph, (gene, _hmm) in clf.models.items()}
    models = []
    for e in evidence:
        if e.role != "core_MAT" or e.status not in MODELLED_STATUSES or not e.exons:
            continue
        span = types.SimpleNamespace(
            contig=e.contig, start=e.start, end=e.end, strand=e.strand, frameshifts=0, cds_blocks=(),
            exons=[types.SimpleNamespace(start=a, end=b) for a, b in e.exons])
        protein = translate(span)
        scores = score_proteins(clf, [protein]) if protein else {}
        models.append(CoreModel(
            gene=e.gene_name, contig=e.contig, start=e.start, end=e.end, protein_aa=len(protein),
            cds_aa=sum(b - a + 1 for a, b in e.exons) // 3,
            scores={class_of[i]: v for i, v in scores.items() if class_of.get(i)}, lost_to=e.lost_to))
    return evaluate(models, max_separation_bp)
