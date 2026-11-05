"""decisions.tsv: one proposed decision per genome of #189 (the 80 genomes with a secondary-system locus in v0.7.0). All decisions are `proposed`: the curator signs off.
Inputs: crosscall_v070.tsv (v0.7.0 flags), evidence.tsv (v0.8.x calls with per-locus evidence), flank10.tsv (flank genes near the loci of the ten two-family genomes).
Tier 1 (policy): v0.8.x reports one family; decision = kept family. Tier 2 (hand): v0.8.x still reports two families; decisions below."""
import csv, collections
cc = [x for x in csv.DictReader(open("crosscall_v070.tsv"), delimiter="\t") if x["suggested_action"].startswith("curate")]
ev = collections.defaultdict(list)
for x in csv.DictReader(open("evidence.tsv"), delimiter="\t"):
    if x["kind"] == "called": ev[x["genome"]].append(x)
fl = collections.defaultdict(dict)
for x in csv.DictReader(open("flank10.tsv"), delimiter="\t"): fl[x["genome"]][x["family"]] = x["flank_hits"]
T2 = {  # genome prefix: (decision, rationale)
 "GCA_003706035.2": ("second_cassette_candidate_fragment", "MATsc alpha locus with SLA2 and SUI1 flanks; MTL A locus (MTLA1, MTLA2) is a 1.8 kb fragment 1.1 kb from the contig start with no flank in range; complementary idiomorph"),
 "GCA_003706045.2": ("second_cassette_candidate_fragment", "as GCA_003706035.2 (same genus): MTL A fragment 1.1 kb from the contig start; MATsc alpha locus has SLA2 and SUI1"),
 "GCA_003708225.3": ("second_cassette_candidate_fragment", "MATsc alpha locus 121-9292 starts 120 bp from the contig start (SLA2, SUI1, DIC1); MTL A locus 638-2122 on another contig, 637 bp from its start, no flank"),
 "GCA_030450195.1": ("second_cassette_candidate_fragment", "MATsc alpha locus with SUI1 and SLA2; MTL A locus 55 kb into another contig with only TGL1 30 kb away (no SLA2 or DIC1): unconfirmed"),
 "GCA_030567515.1": ("second_cassette_candidate_fragment", "MTL A fragment 633-2387 at the start of contig 1; MATsc alpha locus 876 kb on the same contig with SUI1, SLA2, DIC1 (partial locus); complementary idiomorph on one contig, far apart"),
 "GCA_054949465.1": ("second_cassette_supported", "MATsc alpha locus with SUI1, SLA2, DIC1; the MTL A locus 588-2372 has SLA2 (72% identity) 584 bp upstream, i.e. a SLA2-MATa cassette end, a complementary idiomorph; same arrangement as the Ogataea two-cassette genomes"),
 "GCA_003707675.3": ("second_locus_unresolved", "MATsc MATa locus with SLA2, SUI1, DIC1; MTL alpha locus (MTLalpha1/2 at 41-49%) on another long contig, 211 kb from a contig edge, with no SLA2 or DIC1 flank: second cassette or cross-call, needs the alpha genes checked"),
 "GCA_030557395.1": ("second_locus_unresolved", "MATsc MATa locus with SLA2, SUI1, DIC1; MTL alpha locus 4.6 kb from a contig start, no flank in range"),
 "GCA_030579955.1": ("second_cassette_candidate_fragment", "MATsc MATa partial locus with SLA2, SUI1, DIC1; MTL alpha locus 281-1959 at the start of its contig, no flank"),
 "GCA_040365455.1": ("same_idiomorph_two_families", "both loci are alpha (MTL alpha with MTLA2/MTLalpha1/MTLalpha2 and MATsc alpha with MATA1/MATA2/MATALPHA1/MATALPHA2, SLA2 beside it): not complementary; the MATsc locus has the SLA2 flank, the MTL one SUI1 only; probably one locus called by two panels at two places or a second copy"),
}
cols = ["genome", "species", "order", "tier", "decision", "status", "v070_primary", "v070_secondary", "v070_secondary_genes", "v080_loci", "rationale"]
out = []
for x in cc:
    g = x["genome"]; L = ev[g]; fams = sorted({l["family"].split(":")[1] for l in L})
    loci = "; ".join(f"{l['family'].split(':')[1]} {l['idiomorph']} {l['confidence']} {l['contig']}:{l['start']}-{l['end']} [{l['genes']}] id {l['core_identity_min']}-{l['core_identity_max']}{' ' + l['flags'] if l['flags'] else ''}" for l in L)
    if len(fams) == 2:
        key = next(k for k in T2 if g.startswith(k)); dec, why = T2[key]; tier = 2
        why += " | flanks " + " / ".join(f"{k}: {v or 'none'}" for k, v in fl[g].items())
    else:
        tier = 1; dec = f"single_family_{fams[0] if fams else 'none'}"
        why = f"v0.8.x reports one family (cross-family policy and routing); v0.7.0 had {x['primary_family']} plus a second locus of {x['secondary_families']} at another position"
    out.append([g, x["species"], x["order"], tier, dec, "proposed", x["primary_family"], x["secondary_families"], x["secondary_genes"], loci, why])
w = csv.writer(open("decisions.tsv", "w"), delimiter="\t", lineterminator="\n"); w.writerow(cols); w.writerows(out)
print(collections.Counter((r[3], r[4]) for r in out))
