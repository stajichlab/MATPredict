"""Polish at most N admitted clusters per family per genome, ranked before polishing.

Curator's request, 2026-09-26: implement and measure. The study over existing
runs (docs/notes/2026-09-26_mtl-flanks-runtime-and-cluster-limit.md) found that
in the slow panels 2-4% of admitted clusters become calls and wall time is
near-linear in their gene load; a per-family top-6 cap was PROJECTED to halve
median wall time at 0-1 lost calls per panel. Measured on 561 genomes the same
day (docs/notes/2026-09-26_polish-cap-measured-and-serinales-scan.md): cap 6 cut
compute 135.3 -> 74.6 h and lost 4 of 613 calls. Curator's ruling 2026-09-26:
cap 6 is the default; None turns it off.

Rank = what is known before polishing: distinct live genes, then best identity,
then hit count. Per family, because a phylum-fallback genome is searched
against many families at once and one genome-wide top-N starves the rest.
"""
import json

from MATPredict.detect.pipeline import run_pipeline

from tests.detect.test_pipeline import _model, _tblastn, _write_order, _write_record

ORDER = (
    "phylum: P\nloci:\n"
    "  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
    "      - {name: flk1, role: flanking_conserved}\n"
)

# Three clusters of one family, far apart: c1 has 3 genes (best), c2 has 2 at
# high identity, c3 has 2 at low identity.
HITS = [
    ("mfa1", "c1", 100, 200, 80.0), ("pra1", "c1", 300, 400, 80.0), ("flk1", "c1", 500, 600, 80.0),
    ("mfa1", "c2", 100, 200, 70.0), ("pra1", "c2", 300, 400, 70.0),
    ("mfa1", "c3", 100, 200, 40.0), ("pra1", "c3", 300, 400, 40.0),
]


def _run(tmp_path, **kw):
    _write_order(tmp_path, ORDER)
    _write_record(tmp_path)
    polished = []

    def localize(*a, **k):
        return [_tblastn(g, c, s, e, identity=i) for g, c, s, e, i in HITS]

    def polish(*, gene_name, window, **k):
        polished.append((window[0], gene_name))
        start, end = next((s, e) for g, c, s, e, _ in HITS if g == gene_name and c == window[0])
        return _model(gene_name, window[0], start, end)

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize, polish_with_exonerate=polish, polish_with_miniprot=polish,
        evidence_diagnostics_path=tmp_path / "ed.jsonl", **kw,
    )
    return outcome, {c for c, _ in polished}


def test_the_cap_defaults_to_six(tmp_path):
    from MATPredict.detect.pipeline import DEFAULT_MAX_POLISHED_CLUSTERS_PER_FAMILY
    assert DEFAULT_MAX_POLISHED_CLUSTERS_PER_FAMILY == 6
    _, contigs = _run(tmp_path)
    assert contigs == {"c1", "c2", "c3"}    # 3 clusters, under the default cap


def test_none_turns_the_cap_off(tmp_path):
    _, contigs = _run(tmp_path, max_polished_clusters_per_family=None)
    assert contigs == {"c1", "c2", "c3"}


def test_the_cli_defaults_to_six_and_zero_turns_it_off():
    import argparse
    from MATPredict.detect.cli import _polish_cap_from_args, register_subcommands
    parser = argparse.ArgumentParser()
    register_subcommands(parser.add_subparsers(dest="command"))
    default = parser.parse_args(["detect", "--genome", "g", "--out-dir", "o"])
    assert _polish_cap_from_args(default) == 6
    off = parser.parse_args(["detect", "--genome", "g", "--out-dir", "o",
                             "--max-polished-clusters-per-family", "0"])
    assert _polish_cap_from_args(off) is None


def test_a_cap_polishes_only_the_top_ranked_clusters(tmp_path):
    _, contigs = _run(tmp_path, max_polished_clusters_per_family=2)
    assert contigs == {"c1", "c2"}          # c3: fewest genes tied with c2, lower identity


def test_a_capped_cluster_is_not_reported(tmp_path):
    outcome, _ = _run(tmp_path, max_polished_clusters_per_family=1)
    assert {r.contig for r in outcome.results} == {"c1"}


def test_diagnostics_mark_the_capped_clusters(tmp_path):
    _run(tmp_path, max_polished_clusters_per_family=1)
    rows = [json.loads(line) for line in (tmp_path / "ed.jsonl").read_text().splitlines()]
    ev = {r["contig"]: r for r in rows if r["kind"] == "evidence"}
    assert ev["c1"]["admitted"] and not ev["c1"]["polish_capped"]
    assert ev["c2"]["admitted"] and ev["c2"]["polish_capped"]
    assert ev["c3"]["polish_capped"]
