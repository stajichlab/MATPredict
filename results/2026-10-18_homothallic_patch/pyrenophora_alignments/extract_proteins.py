"""Core MAT proteins of chosen genomes from a patched run's reports, for alignment.
python extract_proteins.py RUNS_DIR LIB GENOME_LIST OUT_FAA  (run on the HPCC; GENOME_LIST: group<TAB>genome)"""
import gzip, shutil, sys, tempfile
from pathlib import Path
import yaml
from MATPredict.detect.pipeline import _translate_model

runs, lib, lst, out = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
tmp = Path(tempfile.mkdtemp())
with open(out, "w") as fo:
    for line in open(lst):
        grp, g = line.split()
        d = yaml.load(open(runs / g / "detection_report.yaml"), Loader=L)
        fa = tmp / "g.fna"
        with gzip.open(lib / f"{g}.fa.gz", "rt") as fi, open(fa, "w") as f:
            shutil.copyfileobj(fi, f)
        cache = {}
        for x in d["detected"]:
            if x["family"] != "Ascomycota:MAT":
                continue
            for e in x["gene_evidence"]:
                if e["role"] != "core_MAT" or not e.get("exons"):
                    continue
                span = type("S", (), dict(contig=e["contig"], start=e["start"], end=e["end"], strand=e["strand"], frameshifts=0,
                                          cds_blocks=(), exons=[type("E", (), dict(start=a["start"], end=a["end"])) for a in e["exons"]]))
                p = _translate_model(fa, span, d["genetic_code"] or 1, cache)
                if p:
                    fo.write(f">{grp}|{g}|{e['gene']}|{e['contig']}:{e['start']}-{e['end']}{e['strand']}|{e['status']}|id{e['identity']}|cov{e['coverage']}\n{p}\n")
