// Phase 1 spike: can the MATPredict detect CLI import under Pyodide?
// Usage: node import_check.mjs <path-to-wheel>
import { loadPyodide } from "pyodide";
import { readFileSync } from "node:fs";
import { basename } from "node:path";

const wheel = process.argv[2];
const pyodide = await loadPyodide();
await pyodide.loadPackage(["micropip", "biopython", "pyyaml", "jsonschema", "zstandard", "requests"]);

// micropip will not install a local file path under Node, so unpack the pure-Python
// wheel into a directory on sys.path. Same effect for a py3-none-any wheel.
pyodide.FS.writeFile(`/tmp/${basename(wheel)}`, readFileSync(wheel));
await pyodide.runPythonAsync(`
import sys, zipfile
zipfile.ZipFile("/tmp/${basename(wheel)}").extractall("/site_pkgs")
sys.path.insert(0, "/site_pkgs")
`);

const result = await pyodide.runPythonAsync(`
import importlib, traceback
out = {}
for mod in [
    "MATPredict",
    "MATPredict.detect.search",
    "MATPredict.detect.pipeline",
    "MATPredict.detect.cli",
    "MATPredict.detect.classifier",
    "MATPredict.db.taxonomy",
    "MATPredict.db.local_taxonomy",
]:
    try:
        importlib.import_module(mod)
        out[mod] = "ok"
    except Exception as exc:
        out[mod] = f"{type(exc).__name__}: {exc}"
out
`);
console.log(JSON.stringify(result.toJs({ dict_converter: Object.fromEntries }), null, 2));
