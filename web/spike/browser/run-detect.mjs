// Run `matpredict detect` under Pyodide with every external tool served by WebAssembly modules.
// Node harness standing in for the browser page: the same JS host and Python bridge run in both.
//
// usage: node run-detect.mjs <genome.fasta> <out-dir> [--phylum Ascomycota]
// env:   WHEEL (MATPredict wheel), DB_DIR (curated db/), WASM_DIR (directory with *-mem.js modules),
//        LOCALIZER_DIR (web/spike/localizer), BRIDGE_DIR (web/spike/browser/python)
import { loadPyodide } from "/Users/jstajich/projects/MATPredict/web/spike/node_modules/pyodide/pyodide.mjs";
import { readFileSync, readdirSync, statSync, mkdirSync, writeFileSync } from "node:fs";
import { join, basename, dirname } from "node:path";
import { pathToFileURL } from "node:url";
import { ToolHost } from "./host/tool-host.mjs";
import { runWasiCommand, WORK_DIR } from "./host/wasi-runner.mjs";

const [genome, outDir, ...rest] = process.argv.slice(2);
const phylum = rest[rest.indexOf("--phylum") + 1] ?? "Ascomycota";
const ENGINE = process.env.ENGINE ?? "losat";  // tblastn engine: "losat" or "mmseqs"
const home = "/Users/jstajich/projects/MATPredict";
const WHEEL = process.env.WHEEL ?? join(home, "web/spike/dist/matpredict-0.7.0-py3-none-any.whl");
const DB_DIR = process.env.DB_DIR ?? join(home, "db");
const WASM_DIR = process.env.WASM_DIR ?? join(process.env.HOME, "wasm-build");
const LOC_DIR = process.env.LOCALIZER_DIR ?? join(home, "web/spike/localizer");
const BRIDGE_DIR = process.env.BRIDGE_DIR ?? join(home, "web/spike/browser/python");

const py = await loadPyodide();
await py.loadPackage(["biopython", "pyyaml", "jsonschema", "zstandard", "requests", "micropip"]);

// Mount the database and genome in Pyodide's filesystem (in a browser these come from fetch/upload).
function mount(src, dst) {
  py.FS.mkdirTree(dst);
  for (const name of readdirSync(src)) {
    const p = join(src, name);
    if (statSync(p).isDirectory()) mount(p, `${dst}/${name}`);
    else py.FS.writeFile(`${dst}/${name}`, readFileSync(p));
  }
}
mount(DB_DIR, "/db");
py.FS.mkdirTree("/data");
py.FS.writeFile(`/data/${basename(genome)}`, readFileSync(genome));
py.FS.mkdirTree("/out");

// Package and bridge files into Pyodide.
py.FS.writeFile("/tmp/matpredict.whl", readFileSync(WHEEL));
const PYHMMER_WHL = process.env.PYHMMER_WHL ?? join(home, "web/spike/browser/pyhmmer/pyhmmer-0.12.3-cp312-cp312-pyodide_2024_0_wasm32.whl");
py.FS.writeFile("/" + basename(PYHMMER_WHL), readFileSync(PYHMMER_WHL));
py.FS.mkdirTree("/psutil_stub");
py.FS.writeFile("/psutil_stub/psutil.py", readFileSync(join(home, "web/spike/browser/pyhmmer/psutil.py")));
py.FS.mkdirTree("/browser");
for (const f of ["localizer.py"]) py.FS.writeFile(`/browser/${f}`, readFileSync(join(LOC_DIR, f)));
py.FS.writeFile("/browser/bridge.py", readFileSync(join(BRIDGE_DIR, "bridge.py")));

// Tool factories: each *-mem module is a MODULARIZE build with the in-memory filesystem.
const factories = {
  mmseqs: (await import(pathToFileURL(join(WASM_DIR, "mmseqs-repro/build/src/mmseqs-mem.js")).href)).default,
  exonerate: (await import(pathToFileURL(join(WASM_DIR, "exo-wasm/src/program/exonerate-mem")).href)).default,
  miniprot: (await import(pathToFileURL(join(WASM_DIR, "miniprot-mem2/miniprot.mjs")).href)).default,
};
const LOSAT_WASM = process.env.LOSAT_WASM ?? join(process.env.HOME, "wasm-build/losat/LOSAT/target/wasm32-wasip1/release/LOSAT.wasm");
const wasiTools = { losat: await WebAssembly.compile(readFileSync(LOSAT_WASM)) };
const host = new ToolHost({ pyFS: py.FS, factories, wasiTools, wasiRunner: { runWasiCommand, WORK_DIR } });
globalThis.js_run = (argvJson) => host.run(JSON.parse(argvJson)).then((r) => JSON.stringify(r));

const result = await py.runPythonAsync(`
import os, sys, zipfile, argparse, traceback
os.environ["MATPREDICT_DB_ROOT"] = "/db"
zipfile.ZipFile("/tmp/matpredict.whl").extractall("/site_pkgs")
sys.path.insert(0, "/psutil_stub")
import micropip
await micropip.install("emfs:/${basename(PYHMMER_WHL)}", deps=False)
sys.path[:0] = ["/site_pkgs", "/browser"]
import bridge  # installs the subprocess routing before any MATPredict import
bridge.ENGINE = "${ENGINE}"
from MATPredict.detect import cli as detect_cli
parser = argparse.ArgumentParser(prog="matpredict")
subs = parser.add_subparsers(dest="cmd")
detect_cli.register_subcommands(subs)
args = parser.parse_args(["detect", "--genome", "/data/${basename(genome)}", "--phylum", "${phylum}",
                          "--preset", "v0.7", "--out-dir", "/out"])
error = None
try:
    rc = args.func(args)
except BaseException as exc:
    rc, error = None, f"{type(exc).__name__}: {exc}\\n" + traceback.format_exc(limit=6)
(rc, error)
`);
const [rc, error] = result.toJs();
if (error) console.error("detect failed:", error);
console.log("detect return code:", rc);

// Copy the report back out to the host.
mkdirSync(outDir, { recursive: true });
for (const name of py.FS.readdir("/out")) {
  if (name === "." || name === "..") continue;
  writeFileSync(join(outDir, name), py.FS.readFile(`/out/${name}`));
}
console.log("wrote", outDir);
