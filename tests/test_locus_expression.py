"""Tests for scripts/locus_expression.py."""
import importlib.util
import pathlib

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
spec = importlib.util.spec_from_file_location("locus_expression", SCRIPTS / "locus_expression.py")
le = importlib.util.module_from_spec(spec)
spec.loader.exec_module(le)


def test_parse_bedcov_mean_depth(tmp_path):
    p = tmp_path / "b.tsv"
    p.write_text("c1\t0\t1000\tg1\t5000\nc1\t2000\t2500\tg2\t0\n")
    assert le.parse_bedcov(str(p)) == {"g1": (1000, 5.0), "g2": (500, 0.0)}


def test_percentile_of():
    assert le.percentile_of(5.0, [1.0, 2.0, 5.0, 9.0]) == 0.5
    assert le.percentile_of(0.0, [1.0, 2.0]) == 0.0


def test_locus_transcripts_helpers(tmp_path):
    spec2 = importlib.util.spec_from_file_location("locus_transcripts", SCRIPTS / "locus_transcripts.py")
    lt = importlib.util.module_from_spec(spec2)
    spec2.loader.exec_module(lt)
    gtf = tmp_path / "x.gtf"
    gtf.write_text('c1\tStringTie\ttranscript\t100\t500\t.\t-\t.\ttranscript_id "t1"; gene_id "g1";\n'
                   'c1\tStringTie\texon\t100\t200\t.\t-\t.\ttranscript_id "t1"; gene_id "g1";\n'
                   'c1\tStringTie\texon\t400\t500\t.\t-\t.\ttranscript_id "t1"; gene_id "g1";\n'
                   'c2\tStringTie\texon\t1\t50\t.\t+\t.\ttranscript_id "t2"; gene_id "g2";\n')
    r = lt.read_gtf_transcripts(str(gtf), "c1")
    assert list(r) == ["t1"] and r["t1"]["strand"] == "-" and r["t1"]["exons"] == [(100, 200), (400, 500)]
    assert lt.overlaps(r["t1"]["exons"], 450, 600)
    assert not lt.overlaps(r["t1"]["exons"], 600, 700)
