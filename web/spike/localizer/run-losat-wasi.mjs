import { WASI } from "node:wasi";
import { readFileSync } from "node:fs";
import { openAsBlob } from "node:fs";
const [wasmPath, ...args] = process.argv.slice(2);
const wasi = new WASI({ version: "preview1", args: ["LOSAT", ...args], env: {}, preopens: { "/work": process.cwd() } });
const mod = await WebAssembly.compile(readFileSync(wasmPath));
const inst = await WebAssembly.instantiate(mod, wasi.getImportObject());
const code = wasi.start(inst);
process.exitCode = code;
