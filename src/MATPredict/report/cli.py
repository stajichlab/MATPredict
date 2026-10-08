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
