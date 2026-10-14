"""Route the detect pipeline's external tool calls to WebAssembly, inside Pyodide.

Install before any MATPredict import: the search wrappers bind `subprocess.run` as a default argument
when their module is imported, so replacing it later has no effect.

Routing:
  makeblastdb, tblastn  -> the localizer module (Python), which calls MMseqs2 through the host
  mmseqs                -> host (MMseqs2 module steps)
  miniprot, exonerate   -> host
Anything else is an error: nothing else should be on the detect path.

`js_run` is a JavaScript function set by the host: js_run(json_argv) -> Promise of json
{returncode, stdout, stderr}. run_sync waits on it; this needs JSPI (Node and recent Chromium).
"""
import json
import shlex
import subprocess

from js import js_run  # JavaScript host: js_run(json_argv) -> Promise of json text
from pyodide.ffi import run_sync  # blocks this Python call until the promise settles (JSPI)
import localizer


ENGINE = "losat"  # set by the page or harness: "losat" or "mmseqs"


class CompletedProcess:
    def __init__(self, args, returncode, stdout, stderr):
        self.args = args
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def _host(argv):
    res = json.loads(str(run_sync(js_run(json.dumps(argv)))))
    return CompletedProcess(argv, res["returncode"], res["stdout"], res["stderr"])


def _run_argv(argv, check):
    """Runner for the localizer: raises CalledProcessError like subprocess.run(check=True)."""
    done = _host(argv) if argv[0] in ("mmseqs", "losat") else _unsupported(argv)
    if check and done.returncode != 0:
        raise subprocess.CalledProcessError(done.returncode, argv, output=done.stdout, stderr=done.stderr)
    return done


def _unsupported(argv):
    return CompletedProcess(argv, 127, "", f"not a WebAssembly tool: {argv[0]}")


def run(args, *popenargs, input=None, capture_output=False, text=False, check=False, timeout=None,
        stdout=None, stderr=None, **kwargs):
    argv = list(args) if isinstance(args, (list, tuple)) else shlex.split(args)
    name = argv[0]
    if name == "makeblastdb":
        localizer.makeblastdb(argv[1:])
        done = CompletedProcess(argv, 0, "", "")
    elif name == "tblastn":
        # Engine chosen once per run: LOSAT (default, WASI) or MMseqs2 (translated search).
        if ENGINE == "losat":
            out = localizer.losat_tblastn(argv[1:], lambda a: _run_argv(a, True))
        else:
            out = localizer.tblastn(argv[1:], ["mmseqs"], lambda a: _run_argv(a, True))
        done = CompletedProcess(argv, 0, out, "")
    elif name in ("miniprot", "exonerate", "mmseqs", "losat"):
        done = _host(argv)
    else:
        done = _unsupported(argv)
    if check and done.returncode != 0:
        raise subprocess.CalledProcessError(done.returncode, argv, output=done.stdout, stderr=done.stderr)
    return done


subprocess.run = run
