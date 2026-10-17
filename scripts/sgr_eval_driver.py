"""Evaluate `--strong-single-gene-policy` off vs on on a list of genomes.

For each genome the detect command runs twice (off, then on) in ONE process, with the
search and polish calls memoised across the two runs: they are deterministic functions of
their arguments, so the `on` run is exactly what an independent run would produce, but it
only pays for the extra clusters the rule admits to polishing. The `off` run is never
cached (it runs first).

Usage: sgr_eval_driver.py LIST.tsv OUT_DIR
LIST.tsv columns: genome_id <tab> taxid <tab> fasta_path (plain or .gz)
Writes OUT_DIR/<genome_id>/{off,on}/detection_report.yaml and appends one JSON line per
genome to OUT_DIR/summary.<pid>.jsonl.
"""
import argparse, gzip, inspect, json, os, shutil, sys, tempfile, time, traceback
from pathlib import Path

import yaml

from MATPredict.detect import cli

REAL = cli.run_pipeline
CACHE: dict = {}
SIG = inspect.signature(REAL).parameters


def _memo(name, fn, keyf):
    def w(*a, **k):
        key = (name, keyf(a, k))
        if key not in CACHE:
            CACHE[key] = fn(*a, **k)
        return CACHE[key]
    return w


def _search_key(a, k):
    kk = {x: v for x, v in k.items() if x not in ("runner", "reference_fasta", "record_families")}
    return repr((a[0], [f.key for f in a[1]], sorted(kk.items(), key=lambda t: t[0])))


def _polish_key(a, k):
    return repr((k["gene_name"], k["window"], k["family"].key, k.get("genetic_code")))


def wrapped(*a, **k):
    for name, keyf in (("search_fast_path", _search_key), ("search_localize", _search_key)):
        k.setdefault(name, _memo(name, SIG[name].default, keyf))
    for name in ("polish_with_exonerate", "polish_with_miniprot"):
        k.setdefault(name, _memo(name, SIG[name].default, _polish_key))
    return REAL(*a, **k)


cli.run_pipeline = wrapped


def parser():
    p = argparse.ArgumentParser()
    cli.register_subcommands(p.add_subparsers(dest="cmd"))
    return p


def compact(report_path):
    doc = yaml.safe_load(Path(report_path).read_text())
    loci = []
    for d in doc.get("detected") or []:
        ss = d.get("strong_single_gene")
        loci.append({
            "family": d.get("family"), "contig": d["contig"], "start": d["start"], "end": d["end"],
            "idiomorph": d.get("idiomorph"), "confidence": d.get("confidence"),
            "locus_class": d.get("locus_class"), "genes": [g["gene"] for g in d.get("gene_evidence") or []],
            "strong": ({"identity": ss["core_identity"], "coverage": ss["core_coverage"],
                        "cluster_identity": ss.get("cluster_identity")} if ss else None),
            "ids": [g.get("identity") for g in d.get("gene_evidence") or []],
        })
    return {"loci": loci, "two_idiomorphs": len(doc.get("two_idiomorphs") or []),
            "not_detected": [n.get("reason", "")[:120] for n in doc.get("not_detected") or []][:3],
            "suppressed": len(doc.get("suppressed_loci") or [])}


def main():
    lst, out = Path(sys.argv[1]), Path(sys.argv[2])
    scratch = Path(os.environ.get("SCRATCH") or tempfile.gettempdir())
    P = parser()
    for line in lst.read_text().splitlines():
        if not line.strip():
            continue
        gid, taxid, fasta = line.split("\t")[:3]
        CACHE.clear()
        row = {"genome": gid}
        tmp = Path(tempfile.mkdtemp(dir=scratch, prefix="sgr_"))
        try:
            g = tmp / "genome.fa"
            op = gzip.open if fasta.endswith(".gz") else open
            with op(fasta, "rb") as fi, open(g, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            for pol in ("off", "on"):
                d = out / gid / pol
                t0 = time.time()
                args = P.parse_args(["detect", "--genome", str(g), "--taxid", str(taxid),
                                     "--out-dir", str(d), "--strong-single-gene-policy", pol])
                cli._cmd_detect(args)
                row[pol] = compact(d / "detection_report.yaml")
                row[pol + "_s"] = round(time.time() - t0)
        except Exception:
            row["error"] = traceback.format_exc()[-800:]
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        with open(out / f"summary.{os.getpid()}.jsonl", "a") as fh:
            fh.write(json.dumps(row) + "\n")


if __name__ == "__main__":
    main()
