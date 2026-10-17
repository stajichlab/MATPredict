"""Join the vote notes, the per-gene table, the locus table, the BioSample fetch and the curated verdicts into flip_review.tsv."""
import csv, re
from pathlib import Path
here = Path(__file__).resolve().parent
vote = Path(here.parent / "2026-10-17_partial_coverage_vote")
notes = list(csv.DictReader(open(vote / "vote_notes_on.tsv"), delimiter="\t"))
off = {r["genome"]: r for r in csv.DictReader(open(vote / "campaign_off_on" / "per_genome.tsv"), delimiter="\t")}
for r in csv.DictReader(open(vote / "truth_off_on" / "per_genome.tsv"), delimiter="\t"): off.setdefault(r["genome"], r)
genes, loci = {}, {}
for r in csv.DictReader(open(here / "genes.tsv"), delimiter="\t"): genes.setdefault(r["genome"], {}).setdefault(r["gene"], r)
for r in csv.DictReader(open(here / "loci.tsv"), delimiter="\t"): loci[r["genome"]] = r
bs = {r["assembly"]: r for r in csv.DictReader(open(here / "bs.tsv"), delimiter="\t")}
org = {l.split("\t")[0]: l.strip().split("\t")[1] for l in open(here / "orgs.tsv")}
gw = {}
for r in csv.DictReader(open(here / "gw_boxes.tsv"), delimiter="\t"):
    k = (r["genome"][:15], "alpha" if r["profile"].endswith("1-1") else "hmg")
    gw[k] = max(gw.get(k, 0), float(r["score"]))
truth = {r["accession"]: r for r in csv.DictReader(open(here.parent / "2026-10-17_biosample_truth" / "truth_set.tsv"), delimiter="\t")}
V = {}  # accession -> (verdict, confidence, evidence)
def put(accs, verdict, conf, ev):
    for a in accs.split(): V[a] = (verdict, conf, ev)
BOT_TRUTH_EV = "Botrytis MAT1-2 architecture (full MAT1-2-1 + MAT1-2-10, MAT1-1-1 remnant lacking residues 1-139, no MAT1-1-5); BioSample mating_type 1-2 (submitter label, method not stated) agrees with the flip; same as the 17-genome curation note"
put("GCA_015147945.1 GCA_015147965.1 GCA_015147985.1 GCA_015148035.1 GCA_015148675.1 GCA_015148715.1 GCA_015148775.1 GCA_015148805.1 GCA_015148815.1 GCA_015148845.1 GCA_015148905.1 GCA_015148915.1 GCA_015148935.1 GCA_015149525.1 GCA_015149565.1 GCA_015149655.1 GCA_000292645.1",
    "flip is right", "high", BOT_TRUTH_EV + " (truth set)")
put("GCA_000227075.1 GCA_019186565.1 GCA_040955445.1 GCA_056956565.1 GCA_004786245.1",
    "flip is right", "high", "sampled Botrytis check: same architecture as the truth genomes (MAT1-2-1 100% identical to SAS405 KF944386, MAT1-2-10 present, MAT1-1-1 remnant = reference residues 127/139-353 only, no MAT1-1-5); no BioSample mating statement; GCA_000227075 is the 2011 Genoscope assembly of T4, whose 2013 assembly is a truth MAT1-2 genome (same MAT1-2-1/remnant numbers 100/62 %)")
put("GCA_019186575.1 GCA_032197465.1 GCA_032197565.1 GCA_037039555.1 GCA_052569665.1 GCA_052570965.1 GCA_052570975.1 GCA_052571025.1 GCA_053455415.1",
    "flip is right (not individually examined)", "medium-high", "Botrytis cinerea with the identical gene table (MAT1-2-1 99.5-100%, remnant 62%, no MAT1-1-5, margin 95-98 bits); no BioSample mating statement; only the gene table was inspected, no miniprot ranges")
put("GCA_004379285.1", "flip is right", "medium-high", "Botryotinia calthae: Botrytis MAT1-2 architecture (MAT1-2-1 100%/97% identity to SAS405, MAT1-2-10, MAT1-1-1 remnant = residues 116-353, 1 internal stop, alpha profile 129); SLA2-APN2-COX13 syntenic; no BioSample or literature mating statement")
put("GCA_002162485.1", "flip is right", "medium", "Stromatinia cepivora: MAT1-1-1 remnant (69% coverage, 10 internal stops = degenerate, alpha 131, residues 105-348 only) beside full MAT1-2-1 and MAT1-2-10; no MAT1-1-5. Margin 66 bits. Locus on a 1.6 Mb contig, SLA2/APN2/COX13 syntenic (100 N within 3 kb). No statement of mating type found")
put("GCA_017357885.1", "both idiomorphs genuinely present or fused", "high", "Monilinia vaccinii-corymbosi RL-1: SLA2-MAT1-1-1 (97.7%, alpha 228, intact incl. N-terminus, 9.7 kb intron)-MAT1-2-10-MAT1-2-1 (HMG 188)-MAT1-1-5 (97.6%)-APN2-COX13, one 5.2 Mb chromosome contig; all four MAT genes. Yow et al. 2021 (doi:10.1093/g3journal/jkaa052, same assembly: N50 4.06 Mb) say the MAT locus of the reference genome indicates potential homothallism and a second isolate has a MAT deletion allele. Neither the vote (MAT1-2) nor the flip (MAT1-1) is right; margin 40.5 is between two strong profiles")
put("GCA_029531795.1 GCA_040955825.1 GCA_059058745.1", "both idiomorphs genuinely present or fused", "medium-high", "Sclerotinia nivalis/sativa: MAT1-1-5 intact (100%, 547 bits) + MAT1-2-10 + MAT1-2-1 + MAT1-1-1 pieces in both orientations (N-terminal 1-210 at 96% to S. sclerotiorum, and C-terminal 116-335, the model the rule saw); the S. sclerotiorum Inv+ fused arrangement (Chitrampalam 2013 doi:10.1371/journal.pone.0056895; db record 5180_44ba1). Three independent assemblies (SnTB1, SM8, S. sativa) agree. Flip to MAT1-2 is wrong (the profile margin of 83 bits compares the C-terminal fragment with MAT1-2-1)")
put("GCA_038090695.1", "unclear", "low", "Sclerotinia spermophila CBS 219.46: Botrytis-MAT1-2-like (MAT1-1-1 remnant residues 115-348, alpha 175, 1 stop; MAT1-2-10 and MAT1-2-1 full; no MAT1-1-5 within +-2 kb). But Sclerotinia mating-type switching (S. trifoliorum, Xu 2016 doi:10.1038/srep27083) loses the MAT1-2 genes and leaves MAT1-1, so a MAT1-2-only Sclerotinia would be unexpected; alpha 175 is the highest remnant score; genome-wide search finds no alpha-box ORF above 133 bits and no N-terminal piece (unlike S. nivalis, where an ORF scores 182), which supports the Botrytis reading; no independent statement")
put("GCA_000410735.1 GCA_029298845.1 GCA_029298765.1 GCA_059058445.1 GCA_059058585.1", "flip is right", "medium", "Ophiostoma: full MAT1-2-1 (HMG 160-172, alpha ~2) and a 'MAT1-1-1' that aligns only the C-terminal half (reference residues 315-360 onward) with alpha 5-18 bits (no alpha box). The curated O. montium MAT1-2 record 230073_uamh-11059 carries the same 224-aa MAT1-1-1 remnant beside MAT1-2-1 (alpha 71, aligns 59-224 at 54%). No MAT1-1-2/-3 found. Genome-wide six-frame search: no ORF scores >= 30 bits on the alpha-box profile anywhere in the genome (the HMG box scores 80-87 on the ORFs of the MAT1-2-1 gene). No literature statement for these strains opened" )
put("GCA_001455505.1 GCA_020002355.1 GCA_021396205.1", "flip is right", "medium", "Leptographium lundbergii / Intubia oerlemansii / Sporothrix brunneoviolacea: full MAT1-2-1 (HMG 155-201) and a short 'MAT1-1-1' (17-35% coverage, alpha 0-60). Genome-wide six-frame search: no ORF >= 30 bits on the alpha-box profile in Leptographium and Intubia; Sporothrix has one at the locus only (60 bits, the same short hit). Sporothrix miniprot aligns 660 residues of O. montium MAT1-1-1 at 36% identity, so a divergent MAT1-1-1 homolog cannot be fully excluded. No literature statement opened")
put("GCA_002573605.1", "both idiomorphs genuinely present or fused", "high", "Polytolypa hystricis UAMH7299 (holotype): SLA2-MAT1-2-1 (93.5%, HMG 247)-MAT1-1-1 (alpha 261)-MAT1-1-4-APN2 on one 285 kb contig. Munoz et al. 2018 (doi:10.1038/s41598-018-22816-6, same Broad assembly lineage) report P. hystricis homothallic with both idiomorphs fused at one locus. `undetermined` is the right call for the current labels; `homothallic` would be the correct description")
put("GCA_051941725.1", "both idiomorphs genuinely present or fused", "medium", "Calogaya cf. arnoldii (Kew LIQ161CABI): BioSample says metagenomic assembly from a lichen thallus (metagenomic=1). SLA2-MAT1-2-1 (99.5%, 7 internal stops in the model, HMG 140-324 depending on translation)-MAT1-1-1 (alpha 323)-APN2-COX13 on one 63 kb contig, 5 colinear genes, so a single-genome fused locus is more likely than a mixed assembly; but a metagenome (mycobiont plus other lichen fungi) cannot exclude a chimeric contig. No literature opened for Calogaya (Xanthoria-type homothallism not verified)")
put("GCA_020498755.1", "not a flip (edge-kept)", "low", "Aspergillus fumigatus E-1-75s-2: vote kept MAT1-1 by the edge guard (MAT1-1-1 at 100% identity to A1163 but truncated 30 bp from the end of a 14 kb contig). The same 14 kb contig also carries a complete MAT1-2-1 (100% identical to Af293, HMG 143) with COX13/APN2, i.e. both idiomorph genes on one contig; this looks like a chimeric/collapsed contig or a mixed assembly, not a clean MAT1-1; classified homothallic_candidate by the pipeline. Not examined further")
rows = []
for n in notes:
    g = n["genome"]; a = g[:15]; b = bs.get(a, {}); G = genes[g]; L = loci[g]
    def gv(name, f): return (G.get(name) or {}).get(f, "")
    bsmat = ""
    if b.get("attributes"):
        m = re.search(r"(mating[_ -]?type|mat_type|mating)=([^;]*)", b["attributes"], re.I)
        bsmat = m.group(0) if m else ""
    ver = V[a]
    rows.append(dict(genome=g, organism=org[a], strain=(re.search(r"strain=([^;]*)|isolate=([^;]*)", b.get("attributes", "")) or [None])[0] or "", genus=org[a].split()[0],
        label_off=off[g]["off"], label_on=n["label_on"], decision=n["decision"], winner_gene=n["winner_gene"], winner_cov=n["winner_cov"], other_gene=n["other_gene"], other_cov=n["other_cov"],
        profile_margin=n["profile_margin"], alpha_on_MAT1_1_1=gv("MAT1-1-1", "alpha_best_of_models"), hmg_on_MAT1_2_1=gv("MAT1-2-1", "hmg_best_of_models"),
        MAT1_1_1_internal_stops=gv("MAT1-1-1", "internal_stops"), MAT1_1_5_cov=gv("MAT1-1-5", "coverage"), MAT1_2_10_cov=gv("MAT1-2-10", "coverage"),
        flank_genes=",".join(x for x in ("SLA2", "APN2", "COX13") if x in G), contig_len=L["contig_len"], locus_edge_dist_bp=min(int(L["locus_start"]) - 1, int(L["contig_len"]) - int(L["locus_end"])),
        N_in_3kb=L["n_N_3kb"], biosample=b.get("biosample", ""), biosample_mating=bsmat or ("1-2 (truth_set stated=%s)" % truth[a]["stated"] if a in truth else "none"),
        genome_alpha_box_best_orf=gw.get((a, 'alpha'), 0) if any(k[0] == a for k in gw) else 'not searched', genome_hmg_box_best_orf=gw.get((a, 'hmg'), 0) if any(k[0] == a for k in gw) else 'not searched', verdict=ver[0], confidence=ver[1], evidence=ver[2]))
w = csv.DictWriter(open(here / "flip_review.tsv", "w"), fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
w.writeheader(); w.writerows(rows)
import collections
print(collections.Counter(r["verdict"] for r in rows))
