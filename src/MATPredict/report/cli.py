"""`matpredict report genome`: HTML (and optionally PDF) report of one `detect` run."""
from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from MATPredict import logger
from MATPredict.report.html import render_genome_report

REPORT_YAML = "detection_report.yaml"
_LISTED_OUTPUTS = ("detection_report.yaml", "detected_loci.gff3", "detected_loci.fasta")


def register_subcommands(subparsers) -> None:
    p = subparsers.add_parser("report", help="Readable reports (HTML, PDF) of detect output")
    action = p.add_subparsers(dest="report_action", required=True)
    g = action.add_parser("genome", help="HTML report of one detect run")
    g.add_argument("--run", required=True, help="A `detect --out-dir`, or its detection_report.yaml")
    g.add_argument("--out", default=None, help="HTML file (default: <run>/report.html)")
    g.add_argument("--pdf", default=None, help="Also write a PDF (WeasyPrint or headless Chrome)")
    g.add_argument("--pdf-engine", choices=["weasyprint", "chrome"], default=None,
                   help="PDF engine (default: the first one installed)")
    g.add_argument("--sample", default=None, help="Sample name when the report records none")
    g.set_defaults(func=_cmd_report_genome)
    t = action.add_parser("table", help="genomes.tsv and loci.tsv from a directory of detect runs")
    t.add_argument("--runs", required=True, help="Directory of `detect --out-dir` folders, one per genome")
    t.add_argument("--out-dir", required=True, help="Where genomes.tsv and loci.tsv are written")
    t.add_argument("--samples", default=None, help="CSV/TSV sample sheet; its columns join on the run folder name")
    t.add_argument("--samples-key", default="ASMID", help="Sample-sheet column that matches the run folder names")
    t.set_defaults(func=_cmd_report_table)
    m = action.add_parser("summary", help="Recovery and mating-type counts by group, from the table output")
    m.add_argument("--table-dir", required=True, help="The `report table --out-dir` holding genomes.tsv and loci.tsv")
    m.add_argument("--by", default="phylum", help="Genome-table column to group by (default: phylum)")
    m.set_defaults(func=_cmd_report_summary)


def load_report(path: Path) -> tuple[dict, Path]:
    """The report document and its directory, from a run directory or the YAML itself."""
    path = Path(path)
    yml = path / REPORT_YAML if path.is_dir() else path
    if not yml.is_file():
        raise FileNotFoundError(f"no {REPORT_YAML} at {path}")
    doc = yaml.safe_load(yml.read_text()) or {}
    if not isinstance(doc, dict) or "detected" not in doc:
        raise ValueError(f"{yml} is not a detection report (no `detected` key)")
    return doc, yml.parent


def sample_name(doc: dict, run_dir: Path, given: str | None) -> tuple[str, bool]:
    """(name, came_from_folder): --sample, then the recorded sample, then the genome file name, then the folder."""
    run = doc.get("run") or {}
    genome_file = (run.get("genome") or {}).get("file")
    name = given or run.get("sample") or (genome_file.split(".")[0] if genome_file else None)
    if name:
        return name, False
    return run_dir.resolve().name, True


def write_run_report(run_dir: Path, pdf: bool = False) -> list[Path]:
    """`report.html` (and `report.pdf` when asked) in a finished `detect --out-dir`; the paths written."""
    doc, run_dir = load_report(Path(run_dir))
    sample, from_folder = sample_name(doc, run_dir, None)
    files = [f for f in _LISTED_OUTPUTS if (run_dir / f).exists()]
    kw = {"sample": sample, "sample_from_folder": from_folder, "files": files}
    html_path = run_dir / "report.html"
    html_path.write_text(render_genome_report(doc, **kw), encoding="utf-8")
    written = [html_path]
    if pdf:
        from MATPredict.report.pdf import html_to_pdf
        html_to_pdf(render_genome_report(doc, print_mode=True, **kw), run_dir / "report.pdf")
        written.append(run_dir / "report.pdf")
    return written


def _cmd_report_genome(args: argparse.Namespace) -> int:
    doc, run_dir = load_report(Path(args.run))
    sample, from_folder = sample_name(doc, run_dir, args.sample)
    files = [f for f in _LISTED_OUTPUTS if (run_dir / f).exists()]
    out = Path(args.out) if args.out else run_dir / "report.html"
    kw = {"sample": sample, "sample_from_folder": from_folder, "files": files}
    out.write_text(render_genome_report(doc, **kw), encoding="utf-8")
    print(f"report -> {out}")
    if args.pdf:
        from MATPredict.report.pdf import html_to_pdf
        html = render_genome_report(doc, print_mode=True, **kw)
        engine = html_to_pdf(html, Path(args.pdf), engine=args.pdf_engine)
        logger.info("PDF written with %s", engine)
        print(f"pdf -> {args.pdf} ({engine})")
    return 0


def _cmd_report_table(args: argparse.Namespace) -> int:
    from MATPredict.report.table import build_tables, read_samples, write_tsv
    samples = read_samples(Path(args.samples), args.samples_key) if args.samples else None
    genomes, loci = build_tables(Path(args.runs), samples, args.samples_key)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "genomes.tsv", genomes)
    write_tsv(out / "loci.tsv", loci)
    print(f"table -> {out}/genomes.tsv ({len(genomes)} genomes), {out}/loci.tsv ({len(loci)} loci)")
    return 0


def _cmd_report_summary(args: argparse.Namespace) -> int:
    from MATPredict.report.summary import markdown, read_tsv, summarize
    from MATPredict.report.table import write_tsv
    d = Path(args.table_dir)
    rows = summarize(read_tsv(d / "genomes.tsv"), read_tsv(d / "loci.tsv"), args.by)
    write_tsv(d / f"summary_by_{args.by}.tsv", rows)
    (d / f"summary_by_{args.by}.md").write_text(markdown(rows, args.by), encoding="utf-8")
    print(f"summary -> {d}/summary_by_{args.by}.tsv, {d}/summary_by_{args.by}.md")
    return 0
