# Strong single-gene rule (#97): an independent recheck of the five recovered loci (#261)

Date 2026-11-21. #261 asks for the 12-locus check before the rule is considered. The agent's report of 2026-11-07 (HPCC `results/2026-11-07_issue261/FINDINGS.md`) is on the issue: on 3,100 candidate genomes the rule recovers 5 loci, changes no existing call, one of the five is a defective model. It was not independently re-run. This note re-derives the five loci from the agent's own `recovered.json` (the loci and the rule-on reports) with separate code; it does not repeat the 3,100-genome run, the on/off comparison or the 117-cluster count. Tables: `results/2026-11-21_single_gene_rule_261/` (RUN.txt names the jobs). Nothing changes in any report, call or default; the rule stays off.

## The five loci
| Genome | Species | Gene (idiomorph) | Model | Translated model, Pfam at the gathering cutoff | Flank genes | Contig and locus position |
|---|---|---|---|---|---|---|
| GCA_022559985.1 | *Parastagonospora avenae* | MAT1-1-1 (MAT1-1) | 344 aa, M start, no internal stop, 85.7% identity, 100% coverage | alpha-box PF04769 206.9 bits | none by tblastn | 135 kb contig, locus 3.6 kb from the start |
| GCA_003497985.1 | *P. nodorum* | MAT1-2-1 (MAT1-2) | 345 aa, complete, 100% identity | HMG-box PF00505 68.7, PF09011 25.7 | none | 60.7 kb contig, locus 2.2 kb from the end |
| GCA_003503065.1 | *P. nodorum* | MAT1-2-1 | 345 aa, complete, 100% identity (the same model) | same as the line above | none | 3.6 kb contig, locus 0.3 kb from the start |
| GCA_011034965.1 | *Fusarium oxysporum* | MAT1-2-1 | 223 aa, complete, 91.9% identity | HMG-box PF00505 79.5, PF09011 38.5 | none | 4.5 kb contig, locus 1.5 kb from the start |
| GCA_009802455.1 | *Neurospora discreta* | MAT1-2-1 | 363 aa translated but the CDS is not a multiple of 3, flagged `defective` (frameshift), 93.9% identity, 95% coverage | HMG-box PF00505 63.1, PF09011 33.8 | none | 3.0 kb contig, locus starts at base 2 |

All five are `low`, `idiomorph_gene_only`, and uncalled in the v0.9.0 run without the rule.

## What the recheck confirms
- Four of the five are complete ORFs that carry the domain of their gene at protein level (alpha-box 206.9 bits; HMG-box 63 to 80 bits); the fifth is the defective model, which still carries an HMG-box. These match the agent's table.
- **No flank gene on any of the five contigs.** A tblastn of the 19 SLA2, 24 APN2 and 9 COX13 reference proteins of the db (all Ascomycota records) against each contig, at 1e-5, returns nothing. Three contigs are 3.0 to 4.5 kb, so synteny cannot say anything; the 60.7 kb and 135 kb contigs lack any flank gene entirely. The positive control of the flank search (a called locus finds all three) was not repeated here; the agent reports it.
- **The gate reads no `defective` flag**: the *N. discreta* model passes the rule's 85% identity and 90% coverage gates and has a frameshift. A defective-model gate is the only change the rule needs before use, as the agent said.

## What it adds
- **Two of the five may be one observation**: GCA_003497985.1 and GCA_003503065.1 are two *P. nodorum* assemblies with an identical 345 aa, 100% identity model (two isolates of a clonal species can carry one allele; I did not compare the assemblies). The rule recovers four species, not five.
- **Species**: *P. avenae*, *P. nodorum*, *F. oxysporum*, *N. discreta*; the agent did not tabulate taxonomy.
- Four of the five loci lie within 2.2 kb of a contig end (the *N. discreta* locus starts at base 2) and the fifth is 3.6 kb from one. They look like assembly fragments of a MAT gene whose flanks are missing, which is also why the call has no synteny.

## What it does not settle
- Truth: there is no label for any of the five, so a correct idiomorph call is unverified; the HMG-box is not specific to MAT1-2-1, though the alpha-box (MAT1-1-1) is more so.
- Whether the 117 lone core hits at 85% or more identity that never reach the model gate matter (agent's count, not re-run); the 3,100 genomes cover candidate genomes only.
- The decision to use the rule is the curator's. It stays off.
