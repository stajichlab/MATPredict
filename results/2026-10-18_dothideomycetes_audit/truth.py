"""Truth comparison: BioSample-stated mating types (independent of the pipeline) and strain-matched literature statements.
python truth.py DIR -> tables/truth_comparison.tsv.  The strain-matched rows are NOT independent of the vote for the
reference-record genome (the curated record's own sequence is a reference), but they are independent for the profile built
without that genus, and for the pipeline outcome on the other genomes."""
import sys
from pathlib import Path
import pandas as pd

D = Path(sys.argv[1])
L = pd.read_csv(D / "loci_audit_flagged.tsv", sep="\t", low_memory=False)
G = pd.read_csv(D / "genome_audit_verdict.tsv", sep="\t", low_memory=False)
bs = pd.read_csv(D / "biosample" / "biosample_dothideo_attrs.tsv", sep="\t")
bs = bs[bs.attribute == "mating_type"]
truth = {r.genome: (r.value, "BioSample mating_type " + r.biosample + " strain " + str(r.strain), "independent_of_pipeline_and_db")
         for r in bs.itertuples()}
lit = {
    "GCF_028858645.1_CocheC4_6": ("MAT1-2", "strain C4 = MAT-2 (Turgeon 1993, PMID 8479433)", "strain_matched_literature_reference_in_db"),
    "GCF_000354255.1_CocheC4_1": ("MAT1-2", "strain C4 = MAT-2 (Turgeon 1993, PMID 8479433)", "strain_matched_literature_reference_in_db"),
    "GCA_028891585.1_CocheC5_4m": ("MAT1-1", "strain C5 = MAT-1 (Turgeon 1993, PMID 8479433)", "strain_matched_literature_reference_in_db"),
    "GCA_000338975.1_CocheC5_3": ("MAT1-1", "strain C5 = MAT-1 (Turgeon 1993, PMID 8479433)", "strain_matched_literature_reference_in_db"),
    "GCF_000219625.1_MYCGR_v2.0": ("MAT1-1", "isolate IPO323 = MAT1-1 (Waalwijk 2002, PMID 11929216)", "strain_matched_literature_reference_in_db"),
}
truth.update(lit)
rows = []
for g, (t, src, kind) in truth.items():
    gg = G[G.genome == g]
    ll = L[(L.genome == g) & L.called]
    if gg.empty:
        rows.append(dict(genome=g, truth=t, source=src, kind=kind, status="not_in_run"))
        continue
    r = gg.iloc[0]
    rows.append(dict(genome=g, species=r["species"], truth=t, source=src, kind=kind, status=r["status"],
                     called=r["idiomorphs_called"], primary=r["primary_idiomorph"], confidence=r["primary_confidence"],
                     locus_class=r["primary_class"], profile=r["primary_profile"], domain=r["primary_domain"],
                     verdict=r["verdict"],
                     vote_correct=(r["primary_idiomorph"] == t), profile_correct=(r["primary_profile"] == t),
                     v060=r.get("v060_idiomorphs"), records_run=r.get("records_run_idiomorphs")))
out = pd.DataFrame(rows)
(D / "tables").mkdir(exist_ok=True)
out.to_csv(D / "tables" / "truth_comparison.tsv", sep="\t", index=False)
print(out.to_string())
