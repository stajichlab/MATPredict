import { loadPyodide } from "/Users/jstajich/projects/MATPredict/web/spike/node_modules/pyodide/pyodide.mjs";
import { readFileSync, writeFileSync } from "node:fs";
import { basename } from "node:path";
const py = await loadPyodide();
await py.loadPackage(["micropip"]);
const whl = process.env.WHL;
py.FS.writeFile("/" + basename(whl), readFileSync(whl));
py.FS.writeFile("/sexM.hmm", readFileSync(process.env.HMM));
py.FS.writeFile("/train.faa", readFileSync(process.env.FAA));
py.FS.writeFile("/score.py", readFileSync(process.env.SCORE));
py.FS.mkdirTree("/psutil_stub");
py.FS.writeFile("/psutil_stub/psutil.py", readFileSync(process.env.PSUTIL));
py.runPython("import sys; sys.path.insert(0, '/psutil_stub')");
py.globals.set("WHL_PATH", "emfs:/" + basename(whl));
await py.runPythonAsync(`
import micropip
await micropip.install(WHL_PATH, deps=False)
`);
await py.runPythonAsync(`
import runpy, sys
sys.argv = ["score.py", "/sexM.hmm", "/train.faa", "/scores.txt"]
runpy.run_path("/score.py", run_name="__main__")
`);
writeFileSync(process.env.OUT, py.FS.readFile("/scores.txt", { encoding: "utf8" }));
console.log("wrote", process.env.OUT);
