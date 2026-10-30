#!/bin/bash
# Reference proteins for the S. schoenii gene models: every curated protein of the Saccharomycotina families.
cd "$(dirname "$0")/../.."
cat results/2026-10-29_yeast_routing/pichia_reference_proteins.faa \
    db/Ascomycota/{Pichiales,Saccharomycetales,Dipodascales}/*/proteins.faa > results/2026-10-29_ss_records/refs.faa
grep -c ">" results/2026-10-29_ss_records/refs.faa
