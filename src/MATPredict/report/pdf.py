"""HTML to PDF without a person at a browser.

Engines, in order: WeasyPrint (Python; an optional install) and a headless
Chrome or Chromium on PATH (or the macOS application). The page's print
stylesheet does the layout in both, so the PDF matches "Print, Save as PDF"
in a browser.
"""
from __future__ import annotations

import logging
import shutil
from contextlib import contextmanager
import subprocess
import tempfile
import time
from pathlib import Path

_CHROME_NAMES = ("chromium", "chromium-browser", "google-chrome", "google-chrome-stable", "chrome")
_CHROME_MAC = (
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
)


class NoPdfEngine(RuntimeError):
    pass


@contextmanager
def _quiet_weasyprint():
    """WeasyPrint logs progress (INFO) and expected warnings about CSS it does not support (the screen-only dark
    scheme), from import onwards; that would bury the caller's own output. Errors still show."""
    loggers = [logging.getLogger(n) for n in ("weasyprint", "weasyprint.progress")]
    levels = [lg.level for lg in loggers]
    for lg in loggers:
        lg.setLevel(logging.ERROR)
    try:
        yield
    finally:
        for lg, level in zip(loggers, levels):
            lg.setLevel(level)


def find_chrome() -> str | None:
    for name in _CHROME_NAMES:
        path = shutil.which(name)
        if path:
            return path
    for path in _CHROME_MAC:
        if Path(path).exists():
            return path
    return None


def available_engine(prefer: str | None = None) -> str:
    """'weasyprint' or 'chrome'. `prefer` picks one; raises NoPdfEngine when none is installed."""
    have = {}
    try:
        with _quiet_weasyprint():
            import weasyprint  # noqa: F401
        have["weasyprint"] = True
    except Exception:  # ImportError, or OSError when Pango is missing
        pass
    if find_chrome():
        have["chrome"] = True
    if prefer:
        if prefer not in have:
            raise NoPdfEngine(f"PDF engine {prefer!r} is not available")
        return prefer
    for name in ("weasyprint", "chrome"):
        if name in have:
            return name
    raise NoPdfEngine(
        "no PDF engine found: install WeasyPrint (`pixi add weasyprint` or `pip install weasyprint`) "
        "or put Chrome/Chromium on PATH; or open the HTML report in a browser and use Print, Save as PDF")


def html_to_pdf(html: str, out_pdf: Path, engine: str | None = None, timeout: int = 120) -> str:
    """Write `out_pdf` from the HTML text; returns the engine used."""
    engine = available_engine(engine)
    out_pdf = Path(out_pdf)
    if engine == "weasyprint":
        with _quiet_weasyprint():
            import weasyprint
            weasyprint.HTML(string=html).write_pdf(str(out_pdf))
        return engine
    chrome = find_chrome()
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "report.html"
        src.write_text(html, encoding="utf-8")
        cmd = [chrome, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-pdf-header-footer",
               f"--user-data-dir={tmp}/profile", f"--print-to-pdf={out_pdf.resolve()}", src.resolve().as_uri()]
        if out_pdf.exists():
            out_pdf.unlink()
        # Chrome on macOS can keep running after the PDF is complete, so wait
        # for the file to stop growing rather than for the process to exit.
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        deadline = time.monotonic() + timeout
        last, stable, exited_at = -1, 0, None
        try:
            while time.monotonic() < deadline:
                size = out_pdf.stat().st_size if out_pdf.exists() else -1
                stable = stable + 1 if size > 0 and size == last else 0
                if stable >= 3:
                    break
                # The launcher can exit before the PDF lands (or while it is
                # still being written); give it a few seconds after exit.
                if proc.poll() is not None:
                    exited_at = exited_at or time.monotonic()
                    if time.monotonic() - exited_at > 5:
                        break
                last = size
                time.sleep(0.25)
        finally:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
        if not out_pdf.exists() or out_pdf.stat().st_size == 0:
            err = (proc.stderr.read() if proc.stderr else "") or ""
            raise RuntimeError(f"Chrome PDF render failed ({proc.returncode}): {err.strip()[-400:]}")
    return engine
