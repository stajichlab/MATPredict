"""datasets/reads/sketch_store/baf_family.py: target construction (haplotype collapse, N gaps, label mapping) and the screen-vs-counts summary."""
import csv
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("baf_family", ROOT / "datasets" / "reads" / "sketch_store" / "baf_family.py")
bf = importlib.util.module_from_spec(spec); spec.loader.exec_module(bf)


def write(p, cols, rows):
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def test_targets(tmp_path):
    # bafD = bafA: two strains share haplotype AAAA, one has CCCC; bafB one haplotype; reverse hits and unknown labels ignored
    (tmp_path / "hits.fa").write_text(">bafD|s1|c1\nAAAA\n>bafD|s2|c1\nAAAA\n>bafD|s3|c1\nCCCC\n>bafB|s1|c2\nGGGG\n>bafC|s1|c3\nTTTT\n>bafPseudo|s1|c4\nACAC\n>other|s1|c5\nNNNN\n")
    bf.main(["targets", "--hits", str(tmp_path / "hits.fa"), "--out", str(tmp_path / "t.fa"), "--tsv", str(tmp_path / "h.tsv"), "--gap", "3"])
    rec = dict(bf.read_fasta(tmp_path / "t.fa"))
    assert list(rec) == ["bafA_all", "bafB_all", "bafC_all", "Pseudobaf_all"]
    assert rec["bafA_all"] == "AAAANNNCCCC"                                   # most common haplotype first, gap between
    assert rec["bafB_all"] == "GGGG"
    haps = list(csv.DictReader(open(tmp_path / "h.tsv"), delimiter="\t"))
    assert [(h["haplotype"], h["strains"]) for h in haps][:2] == [("bafA_hap1", "2"), ("bafA_hap2", "1")]


def test_summarize(tmp_path):
    # strains: A1 and A2 have bafD (assembly); A2 only has one old-style hash; B1 has nothing; EX is flagged; Neb strain has bafA+Neb-h1
    calls = [dict(sample_id="local__A1", strain="A1", source="local", excluded="", clone="c1", bafA_call="P", **{"Nebuchadnezzar-h1_call": "P"}),
             dict(sample_id="local__A2", strain="A2", source="local", excluded="", clone="c2", bafA_call="-", **{"Nebuchadnezzar-h1_call": "-"}),
             dict(sample_id="SAMN9", strain="B-1", source="public", excluded="", clone="c3", bafA_call="-", **{"Nebuchadnezzar-h1_call": "-"}),
             dict(sample_id="local__EX", strain="EX", source="local", excluded="qc_mixture", clone="c4", bafA_call="P", **{"Nebuchadnezzar-h1_call": "P"})]
    write(tmp_path / "calls.tsv", ["sample_id", "strain", "source", "excluded", "clone", "bafA_call", "Nebuchadnezzar-h1_call"], calls)
    hits = {"local__A1": (4, 0, 2, 3), "local__A2": (3, 0, 0, 0), "SAMN9": (0, 1, 0, 2), "local__EX": (5, 0, 0, 0)}   # bafA, bafB, bafC, Pseudobaf hits
    rows = [dict(sample_id=s, target=f"{g}_all", target_hashes=13, specific_hashes=12, hits=h, rel_abund=1.1) for s, v in hits.items() for g, h in zip(("bafA", "bafB", "bafC", "Pseudobaf"), v)]
    write(tmp_path / "screen.tsv", ["sample_id", "target", "target_hashes", "specific_hashes", "hits", "rel_abund"], rows)
    write(tmp_path / "counts.tab", ["taxa", "bafB", "bafC", "bafD", "bafPseudo"],
          [dict(taxa="A1", bafB=0, bafC=1, bafD=1, bafPseudo=1), dict(taxa="A2", bafB=0, bafC=0, bafD=1, bafPseudo=0), dict(taxa="B_1", bafB=1, bafC=0, bafD=0, bafPseudo=1),
           dict(taxa="EX", bafB=0, bafC=0, bafD=1, bafPseudo=0), dict(taxa="GONE", bafB=0, bafC=0, bafD=0, bafPseudo=0)])
    (tmp_path / "d").mkdir()
    (tmp_path / "pl.json").write_text(json.dumps(dict(standalone_divergent=dict(genomes=1, bafA_hit=1, hrmA_or_cgnA_in_genome=0))))
    bf.main(["summarize", "--screen", str(tmp_path / "screen.tsv"), "--calls", str(tmp_path / "calls.tsv"), "--counts", str(tmp_path / "counts.tab"),
             "--out", str(tmp_path / "s.json"), "--detail", str(tmp_path / "d"), "--placement", str(tmp_path / "pl.json")])
    s = json.loads((tmp_path / "s.json").read_text()); v = s["validation"]["counts.tab"]
    assert v["taxa"] == 5 and v["matched_taxa"] == 3 and v["pairs"] == 3                      # EX flagged out, GONE not in the store, B_1 matched by normalised name
    assert v["genes"]["bafA"]["2"] == dict(tp=2, fn=0, fp=0, tn=1)                              # A1 (4 hits), A2 (3 hits) both called at >= 2
    assert v["genes"]["bafA"]["4"] == dict(tp=1, fn=1, fp=0, tn=1)                              # stricter cut loses A2
    assert v["genes"]["bafB"]["1"] == dict(tp=1, fn=0, fp=0, tn=2)
    assert v["previous_single_target_bafA"] == dict(tp=1, fn=1, fp=0, tn=1)                     # the old call missed A2
    p = s["population"]
    assert p["kept"] == 3 and p["genes"]["bafA"]["isolates"] == 2 and p["genes"]["bafB"]["isolates"] == 1 and p["genes"]["bafC"]["isolates"] == 1
    assert p["old_vs_new_bafA"] == dict(both=1, new_only=1, old_only=0, neither=1) and p["newly_called"] == dict(isolates=1, neb_h1_present=0)
    assert p["bafA_all_vs_neb_h1"] == dict(both=1, bafA_only=1, neb_only=0, neither=1) and p["by_source"]["local"] == dict(present=2, total=2)
    cl = p["classes"]                                                                           # A1 is cluster-linked (old target), A2 only the all-alleles screen sees
    assert cl["cluster_linked"]["isolates"] == 1 and cl["cluster_linked"]["neb_h1_present"] == 1 and cl["cluster_linked"]["neb_h1_absent"] == 0
    assert cl["standalone_divergent"]["isolates"] == 1 and cl["standalone_divergent"]["neb_h1_present"] == 0 and cl["standalone_divergent"]["neb_h1_absent"] == 1
    assert cl["cluster_linked"]["clone_groups_any_member"] == 1 and cl["cluster_linked"]["largest_clone_group"] == 1 and cl["standalone_divergent"]["isolates_in_groups_over_1"] == 0
    assert p["placement"]["standalone_divergent"]["hrmA_or_cgnA_in_genome"] == 0
    detail = (tmp_path / "d" / "concordance_counts.tsv").read_text()
    assert "A1" in detail and "EX" not in detail                                                # detail names strains; JSON has none
    assert "A1" not in json.dumps(s) and "SAMN9" not in json.dumps(s)
