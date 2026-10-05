# Campaign dashboard and report: guide

What this is: a summary of the whole-clade MATPredict runs on the BFD genomes (code v0.6.0, 2026-10-03 to 05), shared as an interactive dashboard and as a PDF.
Both are generated from tables already committed under `results/`; nothing here is new computation except the Dothideomycete SLA2 test noted below.

## Files
| File | What it is |
|---|---|
| `analysis/2026-10-05_campaign-report.pdf` | The report as a PDF: reading guide, then the four tabs in order (11 pages, Letter). |
| `analysis/2026-10-05_campaign-overview.md` | Markdown overview: campaign table, definitions, call rates, locus size and content, where the reports are, the re-run question. |
| `analysis/2026-10-05_campaign-dashboard-guide.md` | This guide. |
| `analysis/2026-10-05_dothideomycetes-sla2.md` | The Dothideomycete SLA2 test behind the last chart. |
| `results/2026-10-05_campaign_overview/dashboard.html` | The interactive dashboard (an HTML fragment for the artifact viewer; tabs, hover text on bars). |
| `results/2026-10-05_campaign_overview/campaign_report_print.html` | The print version rendered to the PDF. |
| `results/2026-10-05_campaign_overview/make_dashboard.py`, `make_figures.py` | Generators (Python standard library only). Run `python3 make_dashboard.py` to rebuild the HTML; the PDF comes from headless Chrome (`--print-to-pdf`). |
| `results/2026-10-05_campaign_overview/fig1_*.svg` to `fig3_*.svg`, `campaign_summary.tsv` | Standalone charts and the numbers behind the summary tab. |

## Terms
- **Called**: a genome with at least one MAT call. Detection, not ground truth. Most clades lack a curated reference, so uncalled means a reference gap until shown otherwise.
- **Routed by lineage / phylum fallback**: genomes whose order has its own curated record use it; the rest are searched against the whole phylum and call less reliably.
- **Locus class**: full locus (MAT genes plus expected flank structure), partial locus, gene only (the idiomorph gene without a recognised locus), homothallic candidate (both idiomorphs at one locus).
- **PR-only calls**: Basidiomycota genomes whose only call is the pheromone-receptor family; 1,360 come from a motif scan alone (`verification: unverified`).
- **SAC / SCA**: order of SLA2, COX13 and APN2 in Xylariales. SAC keeps the outgroup order with an open interval where a MAT locus sits; SCA has COX13 between SLA2 and APN2.

## Tabs
1. **Summary**: headline counts; call rate per Ascomycota class and Basidiomycota order (marker = rate without PR-only calls); locus class and confidence by campaign.
2. **Exemplars and outliers**: cards with a source file each.
3. **Xylariales synteny**: gene-order maps of seven genomes; genomes versus species (41 of 257 genomes are one species); neighbour adjacency; SLA2 to APN2 distance.
4. **Locus size and content**: Mucoromycotina length and gene content by genus; filamentous Ascomycota length and flank genes; where SLA2 sits relative to the Dothideomycete MAT locus.

## Reading the charts
Bars are drawn to scale from zero. In the range charts the bar is the interquartile range, the dot the median and the line the full range. In the gene-order maps arrows show strand
and the gold box marks an interval of 3.5 kb or more between SLA2 and the next flank gene. Heat tables shade the share of genomes or loci.

## Findings carried by the report
- Call rates are high where a clade has its own record (Ustilaginales, Wallemiales, Agaricales above 98%) and low where it does not (Orbiliomycetes 9%, Dipodascomycetes 29%).
- PR-only calls make fallback Basidiomycota orders look better than they are (Cantharellales 91% to 17%, Trichosporonales 47% to 8% without them).
- Xylariales: two gene-order changes relative to outgroups; the closed MAT interval of the SCA state is shared by 73 species.
- Dothideomycetes: SLA2 is present in every sampled genome but inside the called locus in only 14% (80% in Sordariomycetes); the detachment is lineage-specific.

## Limits
Counts are descriptive and several runs used earlier code. Genome counts reflect uneven strain sampling. Synteny is gene-level. The Dothideomycete test samples 200 of 2,248 loci.
The dashboard uses the Mycotypha and PR-test items as they stood on 2026-10-05; PR #33 was still a draft.
