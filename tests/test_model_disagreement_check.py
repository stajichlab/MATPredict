"""Fixture test of scripts/model_disagreement_check.py (#290, #349): splice motifs, translation, recommendation and the end-to-end run."""
from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path

import yaml

spec = importlib.util.spec_from_file_location("mdc", Path(__file__).resolve().parents[1] / "scripts" / "model_disagreement_check.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

EXON1 = "ATGGCTGCTAAA"                       # M A A K
INTRON = "GT" + "A" * 40 + "AG"              # 44 bp, canonical
EXON2 = "GCTAAAGGTTAA"                       # A K G stop
GENOME = "CCCC" + EXON1 + INTRON + EXON2 + "CCCC"
GOOD = [(5, 16), (61, 72)]                   # exon 1 at 5-16, intron 17-60, exon 2 at 61-72
SHIFTED = [(5, 18), (63, 72)]                # donor moved by two bases: intron 19-62 is noncanonical, still in frame (14 + 10 nt)


def test_translate_and_genetic_codes():
    assert mod.translate("ATGCTGAAA", 1) == "MLK" and mod.translate("ATGCTGAAA", 12) == "MSK"
    try:
        mod.translate("ATG", 4)
    except ValueError:
        pass
    else:
        raise AssertionError("an unsupported code must raise")


def test_a_canonical_model_is_clean():
    m = mod.analyse_model(GENOME, GOOD, "+", 1, ref_aa=7)
    assert m["protein"] == "MAAKAKG" and m["ends_stop"] and m["starts_met"] and m["internal_stops"] == 0 and m["in_frame"]
    assert [i["class"] for i in m["introns"]] == ["GT-AG"] and m["introns"][0]["length"] == 44 and m["ref_ratio"] == 1.0
    assert mod.issues(m) == []


def test_a_shifted_junction_is_noncanonical():
    m = mod.analyse_model(GENOME, SHIFTED, "+")
    assert m["introns"][0]["class"] == "noncanonical" and any("noncanonical" in s for s in mod.issues(m))


def test_minus_strand_reads_the_splice_sites_on_the_gene_strand():
    rc = mod.revcomp(GENOME)
    n = len(GENOME)
    exons = sorted((n - e + 1, n - s + 1) for s, e in GOOD)
    m = mod.analyse_model(rc, exons, "-", 1)
    assert m["protein"] == "MAAKAKG" and m["introns"][0]["class"] == "GT-AG"


def test_internal_stop_frame_and_short_intron_are_issues():
    stop = "ATGTAAAAAGCTTAA"                  # M * K A stop
    m = mod.analyse_model("CC" + stop + "CC", [(3, 17)], "+")
    assert m["internal_stops"] == 1 and any("internal stop" in s for s in mod.issues(m))
    f = mod.analyse_model("CC" + "ATGGCTGC" + "CC", [(3, 10)], "+")
    assert not f["in_frame"] and "length not a multiple of 3" in mod.issues(f)
    short = mod.analyse_model("ATGGCT" + "GT" + "A" * 10 + "AG" + "GCTTAA", [(1, 6), (21, 26)], "+")
    assert short["introns"][0]["short"] and any("under 30" in s for s in mod.issues(short))


def test_recommendation_rules():
    good = mod.analyse_model(GENOME, GOOD, "+")
    bad = mod.analyse_model(GENOME, SHIFTED, "+")
    assert mod.recommend(good, bad)[0] == "primary model"
    assert mod.recommend(bad, good)[0] == "alternate model"
    assert mod.recommend(bad, bad)[0] == "neither clean"
    both = mod.recommend(good, good, ref_aa=7)
    assert both[0] == "both plausible" and "same length" in both[1]
    shared, only_p, only_a = mod.shared_introns(good, bad)
    assert shared == [] and only_p == [(17, 60)] and only_a == [(19, 62)]


def test_run_end_to_end(tmp_path):
    (tmp_path / "genomes").mkdir()
    with gzip.open(tmp_path / "genomes" / "G_S1.fa.gz", "wt") as fh:
        fh.write(">c1 desc\n" + GENOME + "\n")
    (tmp_path / "db" / "Asco" / "Fam" / "r1").mkdir(parents=True)
    (tmp_path / "db" / "Asco" / "Fam" / "r1" / "proteins.faa").write_text(">r1|gene_index=0|name=MAT1-1-3|role=core_MAT\nMAAKAKG\n>r1|gene_index=1|name=OTHER|role=x\nMM\n")
    (tmp_path / "runs" / "S1").mkdir(parents=True)
    entry = {"gene": "MAT1-1-3", "status": "polished_disagree", "contig": "c1", "start": 5, "end": 72, "strand": "+", "method": "exonerate_refine",
             "identity": 99.0, "reference_record": "r1", "exons": [{"start": s, "end": e} for s, e in GOOD],
             "alternate_model": {"contig": "c1", "start": 5, "end": 72, "strand": "+", "method": "miniprot_refine", "identity": 90.0,
                                 "exons": [{"start": s, "end": e} for s, e in SHIFTED]}}
    (tmp_path / "runs" / "S1" / "detection_report.yaml").write_text(yaml.safe_dump(
        {"genetic_code": 1, "detected": [{"family": "Ascomycota:MAT", "confidence": "medium", "gene_evidence": [entry, {"gene": "SLA2", "status": "polished_agree"}]}]}))
    (tmp_path / "genes.tsv").write_text("sample\tgene\tproposed_action\nS1\tMAT1-1-3\tcurator look at the model (tools disagree)\nS1\tOTHER\tno action\n")
    rows = mod.run(tmp_path / "genes.tsv", tmp_path / "runs", tmp_path / "genomes", "G_{sample}.fa.gz", tmp_path / "db", tmp_path / "out")
    assert len(rows) == 1 and rows[0]["recommendation"] == "primary model" and rows[0]["ref_aa"] == 7 and rows[0]["primary_aa"] == 7
    assert rows[0]["only_primary_introns"] == 1 and rows[0]["only_alternate_introns"] == 1
    assert (tmp_path / "out" / "introns.tsv").read_text().count("\n") == 3
    assert (tmp_path / "out" / "models.faa").read_text().startswith(">S1|MAT1-1-3|primary|c1\nMAAKAKG")


def test_reference_conflict_turns_a_motif_recommendation_into_a_conflict():
    good = mod.analyse_model(GENOME, GOOD, "+")          # 7 aa, clean
    bad = mod.analyse_model(GENOME, SHIFTED, "+")        # noncanonical, 7 aa
    # same length: no conflict whatever the identities
    assert mod.recommend(bad, good, 7, 99.0, 80.0)[0] == "alternate model"
    longer = dict(good, aa=12)
    assert mod.recommend(bad, longer, 7, 99.0, 90.0)[0] == "conflict"            # alternate is 9 points lower AND 5 aa farther
    assert mod.recommend(bad, longer, 7, 99.0, 96.0)[0] == "alternate model"     # within the 5-point margin
    assert mod.recommend(bad, longer, 7, None, 90.0)[0] == "alternate model"     # no identity: no conflict
    assert mod.recommend(good, dict(bad, aa=12), 7, 99.0, 90.0)[0] == "primary model"   # chosen = the primary, which is also better: no conflict
    far = dict(good, aa=12)
    assert mod.recommend(far, bad, 7, 80.0, 99.0)[0] == "conflict"               # primary clean but worse; alternate (noncanonical) is better: symmetric
