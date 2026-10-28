# The six wave 1 disagreements with the `detect` call (2026-10-27)

Wave 1 final: 605 runs have a single-allele reads call and a single-allele `detect` call on an assembly of the same BioSample; 599 agree, 6 do not. Each was checked by typing **the assembly itself** with the species panel (`../assembly_panel_check.py`: the assembly's sequence, two copies, goes through `reads-type` as reads), so the assembly is judged by the same k-mer rule as the reads. Output on the HPCC: `results/2026-10-26_v070_submit/asmcheck.out`.

| Run (strain) | Species | Reads call | `detect` on assembly | Panel on the assembly | Reading |
|---|---|---|---|---|---|
| SRR5420328 (CMW 44207, assembly GUA1) | *Dothistroma septosporum* | MAT1-2 (3.0x, breadth 0.50) | MAT1-1 | MAT1-1 (1.3x, flank depth 1.0) | assembly and `detect` agree; the reads do not |
| SRR5420329 (CMW 37965, assembly GRE1) | *D. septosporum* | MAT1-1 (5.8x, breadth 0.66) | MAT1-2 | MAT1-2 (1.0x, breadth 0.49; flank depth 0.7, below the call floor) | assembly and `detect` agree; the reads do not |
| SRR7642325 (14-48) | *Macrophomina phaseolina* | MAT1-1 (6.3x, breadth 0.87, no MAT1-2) | MAT1-2 | MAT1-2 (breadth 0.96) | assembly and `detect` agree; clean reads signal for the other idiomorph |
| SRR31408079 (FPI024) | *Fusarium proliferatum* | MAT1-1 (10x, breadth 0.58, no MAT1-2), a k=21 retyping | MAT1-2 | MAT1-2 (breadth 0.96, no MAT1-1) | assembly and `detect` agree; the reads (k=21 rescue only) do not |
| SRR18739400 (FJII-L3-CM-P1) | *Aureobasidium melanogenum* | MAT1-1 (11.9x, breadth 1.0) | MAT1-2, a partial locus | no panel types the assembly (flank depth below 1 in every lineage) | the `detect` label rests on a partial locus; the reads are the stronger evidence |
| SRR8399124 (UV-GVT) | *Ustilaginoidea virens* | MAT1-1 (2.9x; trace MAT1-2 0.6x) | MAT1-2 | MAT1-1 (breadth 0.53) | panel and reads agree with each other, `detect` is the odd one; the species is homothallic-like (#138) |

## What each means
- **Two runs of one study that look swapped (*Dothistroma*).** Reads of CMW 44207 carry the idiomorph of the CMW 37965 assembly and the reverse; each assembly agrees with its own `detect` call. The two BioSamples and assemblies carry the same strain names as the runs (NCBI assembly summary), so the exchange is inside the submission, either run-to-BioSample or assembly-to-BioSample. The data cannot say which side is wrong; the third piece of evidence would be the published mating types of these isolates (not checked).
- **Real conflicts, not explained (*Macrophomina* 14-48, *F. proliferatum* FPI024).** The assembly and `detect` agree and the reads give a clean, deep signal for the other idiomorph with none of the assembly's. Possible causes: reads and assembly from different isolates or a mix-up, contamination, a collapsed or mis-assembled locus. Not tested: other runs of the same strains (none are in the sample), mapping the reads to the assembly's locus. For FPI024 the call exists only at k=21 (k=31 gave no call), so a k=21 artefact is also possible; `detect` withheld loci of both idiomorphs on one chromosome of that assembly.
- **The `detect` label is the weak side (*Aureobasidium*, *Ustilaginoidea*).** Both are in the homothallic-label weakness already tracked in #138. They should not count against the reads.

## For the site
The species reports list these runs as disagreements. Until each is explained, say so in the report: "disagrees with the assembly call; the assembly checks out against the panel" (4 runs) and "the assembly call rests on a partial locus or an uncertain label" (2 runs). Counting: 599 of 605 stands; two of the six are not evidence against the reads, two look like a swap of the two *Dothistroma* runs, and two are unresolved.
