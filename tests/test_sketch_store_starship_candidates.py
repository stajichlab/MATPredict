"""datasets/reads/sketch_store/starship_candidates.py: tiers, wide table, per-element summary with clone votes."""
import csv
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("starship_candidates", ROOT / "datasets" / "reads" / "sketch_store" / "starship_candidates.py")
sc = importlib.util.module_from_spec(spec); spec.loader.exec_module(sc)


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def write(p, cols, rows):
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def test_tier_of():
    assert sc.tier_of(1.0, "present", "present") == "strong"
    assert sc.tier_of(0.5, "present", "present") == "probable"        # bafA present on one hash only
    assert sc.tier_of(1.0, "present", "partial") == "probable"
    assert sc.tier_of(1.0, "present", "absent") == "discordant"
    assert sc.tier_of(0.0, "absent", "present") == "discordant"
    assert sc.tier_of(0.0, "absent", "partial") is None


def test_tables(tmp_path):
    # s1,s2 = one clone (strong); s3 probable; s4 discordant, flagged; s5 nothing; s6 strong but local
    rows = {"s1": (1.0, "present", 0.9, "present"), "s2": (1.0, "present", 0.8, "present"), "s3": (1.0, "present", 0.4, "partial"),
            "s4": (1.0, "present", 0.0, "absent"), "s5": (0.0, "absent", 0.0, "absent"), "s6": (1.0, "present", 0.7, "present")}
    write(tmp_path / "screen.tsv", ["sample_id", "target", "fraction", "call"],
          [dict(sample_id=s, target=t, fraction=f, call=c) for s, (bf, bc, nf, nc) in rows.items() for t, f, c in (("bafA", bf, bc), ("Nebuchadnezzar-h1", nf, nc))])
    meta = dict(s1=("c1", "public", ""), s2=("c1", "public", ""), s3=("c3", "public", ""), s4=("c4", "public", "qc_mixture"), s5=("c5", "public", ""), s6=("c6", "local", ""))
    write(tmp_path / "samples.tsv", ["sample_id", "source", "study", "country", "strain", "excluded", "clone", "lineage"],
          [dict(sample_id=s, source=src, study="P", country="USA", strain=f"st{s}", excluded=ex, clone=c, lineage="L1") for s, (c, src, ex) in meta.items()])
    sc.main(["--screen", str(tmp_path / "screen.tsv"), "--samples", str(tmp_path / "samples.tsv"), "--out", str(tmp_path)])
    cand = tsv(tmp_path / "bafA_neb_h1_candidates_v2.tsv")
    assert [(r["tier"], r["sample_id"]) for r in cand] == [("strong", "s1"), ("strong", "s2"), ("strong", "s6"), ("probable", "s3"), ("discordant", "s4")]
    assert next(r for r in cand if r["sample_id"] == "s4")["excluded"] == "qc_mixture" and next(r for r in cand if r["sample_id"] == "s6")["source"] == "local"
    wide = tsv(tmp_path / "starship_calls_wide.tsv")
    assert len(wide) == 6 and wide[0]["bafA_call"] == "P" and wide[2]["Nebuchadnezzar-h1_call"] == "p" and wide[4]["bafA_call"] == "-"
    summ = {r["element"]: r for r in tsv(tmp_path / "starship_summary.tsv")}
    assert summ["bafA"]["kept_isolates"] == "5" and summ["bafA"]["present"] == "4" and summ["bafA"]["clone_groups"] == "4"   # s4 flagged; s1+s2 one clone
    assert summ["bafA"]["clones_present"] == "3"                                                                                  # c1, c3, c6
    assert summ["Nebuchadnezzar-h1"]["partial"] == "1" and summ["Nebuchadnezzar-h1"]["clones_present"] == "2"                       # c1, c6 (c3 partial)
