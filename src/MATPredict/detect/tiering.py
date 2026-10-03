"""Per-family confidence tiering -- see spec section
"Boundary calling and confidence tiering" for the rule this encodes. When
any gene in a family is left unpolished (neither miniprot nor exonerate
--refine could confirm it), the tier is capped at Medium regardless of
flanking genes."""
from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, expected_genes_for_idiomorph
from MATPredict.detect.scoring import FamilyScore

_TIER_DOWNGRADE = {"high": "medium", "medium": "low", "low": "low"}

#: Curator's ruling, 2026-09-27 (variant B' with the closeness guard; replay
#: in results/2026-09-27_tier_rule_replay). A gene whose
#: `present_in_idiomorphs` excludes the called idiomorph is ignored by the
#: tier only when it is weaker than this ...
ALLELE_ABSENT_MAX_IDENTITY = 50.0
#: ... AND at least this many identity points below the called allele's best
#: MODELLED core gene. Without the 50% cap the replay promoted 108 calls on a
#: strong "absent" gene (S. cerevisiae MATA2 at ~100% by a2/alpha2 homology,
#: Rhizopus btbA at 66-99%, collapsed a/alpha Serinales loci); without the gap
#: it promoted calls whose ignored gene was as strong as, or stronger than,
#: the called one (under 10 points apart in 78 of 167 B' risers).
ALLELE_ABSENT_MIN_GAP = 10.0


def allele_absent_genes_to_ignore(
    family: Family,
    idiomorph: str | None,
    identities: dict[str, float | None],
    modelled: set[str] | frozenset[str],
    locus_class: str | None = None,
) -> frozenset[str]:
    """Genes the tier should not count against a call of `idiomorph`.

    A cross-hit to the OTHER allele's gene is expected at a real
    single-allele locus: MAT1-1 and MAT1-2 genes, or the two Wallemia
    receptors, share domains. Counting it did two wrong things. An unmodelled
    cross-hit capped the call at medium through `any_gene_unpolished`, and a
    found cross-hit made the found genes name both alleles, so
    `expected_genes_for_idiomorph` fell back to the full roster and a gene
    the called allele never carries (Wallemia v2's SXI1) counted as missing.
    All 18 Wallemia v2 calls were medium for the second reason alone.

    Only a WEAK, CLEARLY SEPARATED cross-hit is ignored (see the two
    constants). Undetermined or multi-allele calls, calls with no modelled
    core gene of the called allele, and genes with no recorded identity
    ignore nothing -- the guard has nothing to measure against, so the
    existing behaviour stands. A gene with no `present_in_idiomorphs`
    belongs to every allele and is never ignored.
    """
    if not idiomorph or idiomorph == "undetermined":
        return frozenset()
    # Curator's ruling 2026-09-27: a `homothallic_candidate` rests on BOTH
    # alleles by definition, so neither may be dropped from its tier. Case
    # that forced it: Serinales GCA_030462985.1 was raised to high by
    # ignoring its MTLA2 (38.75%) (results/2026-09-27_tier_rule_implemented).
    if locus_class == "homothallic_candidate":
        return frozenset()
    pin = {g["name"]: set(g.get("present_in_idiomorphs") or ()) for g in family.genes}
    if not any(idiomorph in allowed for allowed in pin.values()):
        return frozenset()  # e.g. "a+alpha", or a pattern allele number
    called_best = max(
        (
            float(identities[g["name"]])
            for g in family.genes
            if g["role"] == "core_MAT"
            and g["name"] in modelled
            and identities.get(g["name"]) is not None
            and (not pin[g["name"]] or idiomorph in pin[g["name"]])
        ),
        default=None,
    )
    if called_best is None:
        return frozenset()
    absent = {
        name: identity
        for name, identity in identities.items()
        if pin.get(name) and idiomorph not in pin[name]
    }
    # The guard is applied to the CALL, not gene by gene: if any other-allele
    # gene here is strong, close to the called allele, or of unknown identity,
    # nothing is ignored. Applied per gene, the real-code replay promoted 15
    # calls that kept a strong other-allele gene -- 6 Serinales alpha calls
    # still carrying a modelled MTLA1 (45%) or MTLA2 (58%), a Periconia
    # MAT1-2 call whose kept MAT1-1-3 (39.5%) beat its own best core gene
    # (31.8%) -- exactly the collapsed-locus shape the guard exists to keep at
    # medium (results/2026-09-27_tier_rule_implemented).
    if any(
        identity is None
        or float(identity) >= ALLELE_ABSENT_MAX_IDENTITY
        or called_best - float(identity) < ALLELE_ABSENT_MIN_GAP
        for identity in absent.values()
    ):
        return frozenset()
    return frozenset(absent)


def cap_at_medium(tier: str) -> str:
    """Lower `tier` to `medium` if it is higher, never raise it.

    Distinct from `_TIER_DOWNGRADE`, which moves every tier down one step --
    that would send `medium` to `low`, and `low` means something specific here
    ("an isolated single hit"). A `partial_locus` carrying three genes and a
    conserved flank is genuinely more than an isolated hit, so it must not be
    collapsed into the same bucket.

    Used for `partial_locus`, per the curator's ruling of 2026-09-21: `high` is
    the report's assertion that something IS a MAT locus, and a call that
    cleared the admission floor only by tying it, or that came from the relaxed
    pass, must not make that assertion.
    """
    return "medium" if tier == "high" else tier


def has_flanking_conserved(family: Family) -> bool:
    return any(g["role"] == "flanking_conserved" for g in family.genes)


def _searchable_flanking_conserved(family: Family, score: FamilyScore) -> bool:
    """Does this family have a `flanking_conserved` gene the run could actually
    have found?

    The flanking branch below caps a call at Medium when a declared flank is
    absent. That is only evidence of absence if the run looked for it. A flank
    with no reference protein, or one the curator has taken off the search list
    (`exclude_from_search`), lands in `genes_not_searchable` and must not cap
    anything -- the same reasoning `assign_tier` already applies to the core
    requirement, and the same bug shape as counting `genes_not_searchable`
    against a locus.

    An `optional` flank is excused for a different reason: it is a BONUS, not
    a requirement. Curator's ruling, 2026-09-21 -- "perhaps sla2 is a bonus
    flank to search for but not to penalize if not present". `sla2` is the
    case that needs it: it is genuinely adjacent to MAT in Kluyveromyces
    lactis and Lachancea thermotolerans (sla2 -> MATA1 -> MATA2), and is not
    on chromosome III at all in Saccharomyces -- and one `MATsc` family spans
    both. Marked optional it is searched everywhere, corroborates a call to
    High where the architecture has it, and costs nothing where it does not.
    Finding it still promotes the tier; only its ABSENCE stops mattering.

    Load-bearing for Saccharomyces. The curator ruled on 2026-09-21 that this
    clade has no usable flanking GENE: S. cerevisiae's three cassettes
    (HML/MAT/HMR) are told apart by flanking DNA -- the X/Y/Z homology boxes --
    not by neighbours, and the two flanks MATsc previously declared were
    measured wrong on the S288C reference (SLA2 is not on chromosome III at
    all; CHA1 abuts HML, not MAT). Both are now excluded from search, and
    without this a perfect two-gene MATsc call would be capped at Medium for
    missing a flank nothing ever searched for.
    """
    unsearchable = set(score.genes_not_searchable)
    return any(
        g["role"] == "flanking_conserved"
        and g["name"] not in unsearchable
        and not g.get("optional")
        for g in family.genes
    )


def assign_tier(
    score: FamilyScore,
    family: Family,
    cluster: GeneCluster,
    any_gene_unpolished: bool,
    fragmented: bool,
    ignore_genes: frozenset[str] = frozenset(),
) -> str:
    # `ignore_genes` (from `allele_absent_genes_to_ignore`) leave the found
    # set before the expected core is derived; the caller also drops them
    # from `any_gene_unpolished`. Empty by default: unchanged behaviour.
    # A gene the run held no reference protein for is dropped from the core
    # requirement, exactly as `scoring.score_cluster` drops it from the
    # `fraction_found` denominator. Counting it here would state that the
    # genome failed to show a gene nothing could ever have found -- evidence
    # of absence manufactured from a gap in the reference database.
    #
    # Measured on Schizophyllum commune H4-8 (GCF_000143185.2, genome-only,
    # 2026-09-20): Bbeta was recovered completely -- all 8 curated genes,
    # `genes_missing=[]`, `fraction_found=1.0` -- and was still capped at
    # Medium, because the Bbeta roster's `pheromone_receptor` alias has no
    # reference protein anywhere in `db/`. Balpha (3/3 found) was demoted the
    # same way. Those were the two most complete calls in that run.
    # An `optional: true` gene is excused from the core requirement for the
    # same reason `genes_not_searchable` is, one clause below: requiring it
    # states that the genome failed to show something that was never required.
    # `scoring.score_cluster` already drops optional genes from
    # `fraction_found`; before 2026-09-21 this function did not, so marking a
    # gene optional removed it from the score while still letting its absence
    # cap the tier.
    #
    # Load-bearing for the Cryptococcus recuration: the curator ruled the
    # recombination-trapped genes (STE20, RPO41, RPL39) and the homeodomain
    # genes (SXI1, SXI2) in as optional, because they are expected in
    # Cryptococcus but not established as universal across Tremellales. The
    # required core is the pheromone precursors and the receptor.
    expected_core = {
        g["name"]
        for g in expected_genes_for_idiomorph(
            family, [n for n in score.genes_found if n not in ignore_genes]
        )
        if g["role"] == "core_MAT" and not g.get("optional")
    }
    core_genes = expected_core - set(score.genes_not_searchable)
    core_requirement_relaxed = core_genes != expected_core
    core_found = core_genes.issubset(set(score.genes_found))

    # Spec: Low is "a single gene hit with no other expected genes from the
    # same family found nearby". Keying Low solely on fraction_found == 0
    # made the tier unreachable in practice, because score_cluster never
    # emits a FamilyScore for a family with zero hits -- so a lone isolated
    # hit was indistinguishable from a substantial partial match. An
    # isolated single hit in a family that expects more than one gene is
    # therefore Low; any richer partial match stays Medium.
    isolated_single_hit = len(score.genes_found) <= 1 and len(family.genes) > 1

    if not core_found:
        tier = "low" if score.fraction_found == 0 or isolated_single_hit else "medium"
    elif isolated_single_hit and core_requirement_relaxed:
        # `core_found` can now be satisfied by ONE gene, when every other core
        # gene of the family is unsearchable. Without this branch that single
        # hit would go straight to High in any family with no
        # `flanking_conserved` gene -- which is every Basidiomycota family.
        # One gene is not a locus, however complete the searchable roster
        # technically was.
        #
        # Guarded on `core_requirement_relaxed` so this only ever fires where
        # the relaxation above put it: a lone core gene in a family whose
        # roster was NOT relaxed keeps its existing tier, which for a flanked
        # family with its flank missing is Medium.
        tier = "low"
    elif any_gene_unpolished:
        tier = "medium"
    elif _searchable_flanking_conserved(family, score):
        flanking_found = any(
            g["role"] == "flanking_conserved" and g["name"] in score.genes_found
            for g in family.genes
        )
        tier = "high" if flanking_found else "medium"
    else:
        tier = "high"

    if fragmented and tier != "low":
        tier = _TIER_DOWNGRADE[tier]
    return tier
