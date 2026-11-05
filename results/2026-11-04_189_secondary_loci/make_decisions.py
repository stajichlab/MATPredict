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
 "GCA_003706035.2": ("second_cassette_candidate_fragment", "MATsc alpha locus with SLA2 and SUI1; MTL A locus (MTLA1 is MATa1, MTLA2 the HMG MATa2) is a 1.8 kb fragment 1.1 kb from the contig start; the contig shares no sequence with the alpha contig and overlaps another contig by 199 bp, so it is not a dead end and not an allele (Opus review 2026-11-04: plausible second cassette, unconfirmed; there is no second SLA2 in the genome); next check: what contig PPKF02000034 carries and whether reads span 33 to 34"),
 "GCA_003706045.2": ("second_cassette_candidate_fragment", "as GCA_003706035.2 (same genus): MTL A fragment 1.1 kb from the contig start; MATsc alpha locus has SLA2 and SUI1"),
 "GCA_003708225.3": ("second_cassette_candidate_fragment", "MATsc alpha locus 121-9292 starts 120 bp from the contig start (SLA2, SUI1, DIC1); MTL A locus 638-2122 on another contig, 637 bp from its start, no flank"),
 "GCA_030450195.1": ("second_cassette_candidate_fragment", "MATsc alpha locus with SUI1 and SLA2; MTL A locus 55 kb into another contig with only TGL1 30 kb away (no SLA2 or DIC1): unconfirmed"),
 "GCA_030567515.1": ("second_cassette_candidate_fragment", "MTL A fragment 633-2387 at the start of contig 1; MATsc alpha locus 876 kb on the same contig with SUI1, SLA2, DIC1 (partial locus); complementary idiomorph on one contig, far apart"),
 "GCA_054949465.1": ("allelic_pair", "the two loci are one heterozygous MAT locus, not two cassettes: the 3.5 kb MTL A contig (JBSXOT010000274.1; 3,507 bp, the A locus itself is 2.4 kb) has SLA2 at one end and DIC1 sequence at the other, and both ends match the MATsc alpha contig (JBSXOT010000020.1) flanking its alpha locus at 97.2% (533 bp, 41,292-41,824) and 97.3% (147 bp, 45,203-45,347): an a and an alpha idiomorph between the same SLA2 and DIC1 (Opus review 2026-11-04, BLAST rows confirmed). Weakest link: only 680 bp of shared flank; reads mapped to each idiomorph would confirm half coverage each"),
 "GCA_003707675.3": ("second_cassette_supported", "MATsc MATa locus (DIC1, a2, a1, SLA2, SUI1 on contig 4, mid-contig) and an MTL alpha locus on contig 12, 211 kb from an edge (alpha2 and alpha1 on one strand 5 kb apart, plus a 35 aa SLA2 3-prime fragment 2 kb beyond alpha1 that the 50 kb flank search missed); the two regions share no nucleotide sequence, so this is not an allele pair; independent review (Opus, 2026-11-04) calls it the strongest second-cassette case. Weakest link: which cassette is expressed is unknown; a sister-species or long-read check would confirm"),
 "GCA_030557395.1": ("second_locus_unresolved", "MATsc MATa locus with SLA2, SUI1, DIC1; MTL alpha locus 4.6 kb from a contig start, no flank in range"),
 "GCA_030579955.1": ("second_cassette_candidate_fragment", "MATsc MATa partial locus with SLA2, SUI1, DIC1; MTL alpha locus 281-1959 at the start of its contig, no flank"),
 "GCA_040365455.1": ("duplicate_haplotype", "one locus on two haplotypes, not a cassette pair: contig BCIU01000086.1 is collinear with BCIU01000040.1 along its whole length (92-98%, the alpha1-alpha2 region 98.1%, identical SUI1); the assembly is 29.4 Mb (about 19 Mb expected) with about 46 divergent haplotype-copy contigs, so a partly separated hybrid or diploid. Both copies carry an HMG plus alpha1 and alpha2; the a1 hit (e-7) and a SLA2 truncated at aa 501 occur only on contig 40. Independent review (Opus, 2026-11-04); weakest link is the gene content of each haplotype (translate the a1 region on contig 40, map reads)"),
}
NOTE = "gene names differ by family: MTLA1 is MATa1, MTLA2 is the HMG MATa2 (Pichia/Komagataella type), MTLalpha1/2 are MATalpha1/2; the MATsc MATA2 label can land on the alpha2 homeodomain"
import os
FP = collections.defaultdict(list)
if os.path.exists("flank_pairs.tsv"):
    for x in csv.DictReader(open("flank_pairs.tsv"), delimiter="\t"): FP[x["genome"]].append(f"{x['status']}:{x['aligned_bp']}bp@{x['best_identity']}" + (" same-family" if x["same_family"] == "True" else ""))
cols = ["genome", "species", "order", "tier", "decision", "status", "flank_check", "v070_primary", "v070_secondary", "v070_secondary_genes", "v080_loci", "rationale"]
out = []
for x in cc:
    g = x["genome"]; L = ev[g]; fams = sorted({l["family"].split(":")[1] for l in L})
    loci = "; ".join(f"{l['family'].split(':')[1]} {l['idiomorph']} {l['confidence']} {l['contig']}:{l['start']}-{l['end']} [{l['genes']}] id {l['core_identity_min']}-{l['core_identity_max']}{' ' + l['flags'] if l['flags'] else ''}" for l in L)
    if len(fams) == 2:
        key = next(k for k in T2 if g.startswith(k)); dec, why = T2[key]; tier = 2
        why += " | flanks " + " / ".join(f"{k}: {v or 'none'}" for k, v in fl[g].items())
    else:
        tier = 1; dec = f"kept_family_{fams[0] if fams else 'none'}_by_policy"
        why = f"SUPPRESSED BY POLICY, not independently resolved: v0.8.x reports one family (cross-family policy and routing); v0.7.0 had {x['primary_family']} plus a second locus of {x['secondary_families']} at another position"
    out.append([g, x["species"], x["order"], tier, dec, "proposed", "; ".join(FP.get(g, [])) or "no complementary pair", x["primary_family"], x["secondary_families"], x["secondary_genes"], loci, why])
w = csv.writer(open("decisions.tsv", "w"), delimiter="\t", lineterminator="\n"); w.writerow(cols); w.writerows(out)
print(collections.Counter((r[3], r[4]) for r in out))
