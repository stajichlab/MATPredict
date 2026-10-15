// Phase 1 spike, part 2: where does `matpredict detect` first call an external tool?
// The subprocess.run seam is replaced before any MATPredict import (the wrappers bind
// it as a default argument at import time). Every call is logged and returns a failure,
// so the run stops at the first real tool call.
import { loadPyodide } from "pyodide";
import { readFileSync, readdirSync, statSync } from "node:fs";
import { join, relative, basename } from "node:path";

const REPO = new URL("../../", import.meta.url).pathname;
const wheel = join(REPO, "web/spike/dist/matpredict-0.6.1-py3-none-any.whl");
const dbDir = join(REPO, "db");
const genome = join(REPO, "testset/MAT_loci3.fasta");

const pyodide = await loadPyodide();
await pyodide.loadPackage(["micropip", "biopython", "pyyaml", "jsonschema", "zstandard", "requests"]);

// Mount db/ and the input genome into Pyodide's filesystem.
function mountTree(src, dst) {
  pyodide.FS.mkdirTree(dst);
  for (const name of readdirSync(src)) {
    const p = join(src, name);
    const target = `${dst}/${name}`;
    if (statSync(p).isDirectory()) mountTree(p, target);
    else pyodide.FS.writeFile(target, readFileSync(p));
  }
}
mountTree(dbDir, "/db");
pyodide.FS.mkdirTree("/data");
pyodide.FS.writeFile("/data/" + basename(genome), readFileSync(genome));
pyodide.FS.writeFile("/tmp/" + basename(wheel), readFileSync(wheel));

const out = await pyodide.runPythonAsync(`
import os, sys, zipfile, json, subprocess, argparse, traceback
os.environ["MATPREDICT_DB_ROOT"] = "/db"
zipfile.ZipFile("/tmp/${basename(wheel)}").extractall("/site_pkgs")
sys.path.insert(0, "/site_pkgs")

calls = []
class _Stub:
    def __init__(self, argv, returncode):
        self.args, self.returncode, self.stdout, self.stderr = argv, returncode, "", "stub: not installed"
def _logged_run(cmd, *a, **kw):
    argv = list(cmd) if isinstance(cmd, (list, tuple)) else [cmd]
    calls.append(argv)
    # Refuse the call: the probe only needs to see the first tool invocation.
    return _Stub(argv, 127)
subprocess.run = _logged_run   # must precede the MATPredict imports below

from MATPredict.detect import cli as detect_cli

parser = argparse.ArgumentParser(prog="matpredict")
subs = parser.add_subparsers(dest="cmd")
detect_cli.register_subcommands(subs)
args = parser.parse_args([
    "detect", "--genome", "/data/MAT_loci3.fasta", "--phylum", "Ascomycota",
    "--out-dir", "/data/out",
])

error = None
try:
    rc = args.func(args)
except BaseException as exc:
    rc = None
    error = f"{type(exc).__name__}: {exc}"
    error_tb = traceback.format_exc(limit=6)

json.dumps({"rc": rc, "error": error, "first_tool_calls": calls[:5],
            "tb": error_tb if error else None})
`);
console.log(out);
