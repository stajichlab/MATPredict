"""The pre-publish site check: it passes a complete site and names each kind of break."""
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("check_site", REPO / "web/spike/page/check_site.py")
cs = importlib.util.module_from_spec(spec)
sys.modules["check_site"] = cs
spec.loader.exec_module(cs)


def _site(tmp_path: Path, app: str = 'import { a } from "./a.mjs";\nawait fetch("../browser/x.py");\n') -> Path:
    s = tmp_path / "site"
    files = {
        "index.html": "x", "results.html": "x", "about.html": "x", "data/summary.json": "{}",
        "page/index.html": "x", "page/app.mjs": app, "page/a.mjs": "export const a = 1;",
        "browser/x.py": "x", "page/_assets/matpredict.whl": "w", "page/_assets/bundle/pyhmmer-1.whl": "w",
        # the Guess button: its modules and GmG's WebAssembly (the .wasm is loaded by the .js, no module names it)
        "page/gmg_guess.mjs": "export const g = 1;", "page/gmg_route.mjs": "export const r = 1;", "page/gmg_worker.mjs": "export const w = 1;",
        "page/_assets/gmg/gmg_wasm.js": "export default 1;", "page/_assets/gmg/gmg_wasm_bg.wasm": "w",
    }
    blob = b"archive bytes"
    (s / "page/_assets").mkdir(parents=True, exist_ok=True)
    for rel, text in files.items():
        p = s / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    (s / "page/_assets/db.tar.gz").write_bytes(blob)
    (s / "page/_assets/db-manifest.json").write_text(json.dumps({
        "db": ["x"], "pyhmmerWheel": "pyhmmer-1.whl", "dbArchive": "db.tar.gz",
        "dbArchiveBytes": len(blob), "dbArchiveSha256": hashlib.sha256(blob).hexdigest()}))
    return s


def test_a_complete_site_has_no_problems(tmp_path):
    assert cs.check(str(_site(tmp_path))) == []


def test_a_page_module_that_was_not_copied_is_named(tmp_path):
    s = _site(tmp_path)
    (s / "page/a.mjs").unlink()
    assert cs.check(str(s)) == ["page/app.mjs refers to ./a.mjs, which is not in the site"]


def test_a_fetched_file_that_is_missing_is_named(tmp_path):
    s = _site(tmp_path)
    (s / "browser/x.py").unlink()
    assert cs.check(str(s)) == ["page/app.mjs refers to ../browser/x.py, which is not in the site"]


def test_a_missing_required_file_is_named(tmp_path):
    s = _site(tmp_path)
    (s / "page/_assets/matpredict.whl").unlink()
    assert "missing page/_assets/matpredict.whl" in cs.check(str(s))


def test_a_missing_gmg_wasm_is_named(tmp_path):
    s = _site(tmp_path)
    (s / "page/_assets/gmg/gmg_wasm_bg.wasm").unlink()
    assert "missing page/_assets/gmg/gmg_wasm_bg.wasm" in cs.check(str(s))


def test_an_archive_that_does_not_match_its_manifest_is_named(tmp_path):
    s = _site(tmp_path)
    (s / "page/_assets/db.tar.gz").write_bytes(b"archive bytez")
    problems = cs.check(str(s))
    assert any("does not match the SHA-256" in p for p in problems)
    assert not any("bytes, the manifest says" in p for p in problems)  # same length, only the checksum differs


def test_a_missing_archive_or_pyhmmer_wheel_is_named(tmp_path):
    s = _site(tmp_path)
    (s / "page/_assets/bundle/pyhmmer-1.whl").unlink()
    problems = cs.check(str(s))
    assert "the manifest names the pyhmmer wheel pyhmmer-1.whl, which is not in the bundle" in problems


def test_a_script_an_html_page_imports_that_was_not_copied_is_named(tmp_path):
    s = _site(tmp_path)
    (s / "reads").mkdir()
    (s / "charts.js").write_text("x")
    (s / "reads/species.html").write_text('<link rel="stylesheet" href="../site.css"><script type="module">import { a } from "../charts.js"; import { b } from "../sortable.js";</script>')
    (s / "site.css").write_text("x")
    assert cs.check(str(s)) == ["reads/species.html refers to ../sortable.js, which is not in the site"]


def test_template_strings_and_urls_are_not_treated_as_files(tmp_path):
    s = _site(tmp_path, app='const u = `${ASSETS}/db.tar.gz`; import x from "https://cdn.example/x.mjs"; fetch("examples/e.json");\n')
    assert cs.check(str(s)) == []


def test_main_exit_codes(tmp_path, capsys):
    s = _site(tmp_path)
    assert cs.main(["check_site.py", str(s)]) == 0
    assert "site check: ok" in capsys.readouterr().out
    (s / "page/a.mjs").unlink()
    assert cs.main(["check_site.py", str(s)]) == 1
    assert "which is not in the site" in capsys.readouterr().out


def test_every_file_the_html_pages_refer_to_exists_in_the_repository():
    """The landing, results, About and reads pages, as checked in: each stylesheet, script and import is a file
    in web/intro, apart from generated data under data/ and the browser tool under page/ (assembled from
    web/spike by build-site.sh)."""
    intro = REPO / "web/intro"
    missing = [f"{name}: {ref}" for name, ref in cs.html_references(str(intro))
               if "/data/" not in ref and "/page/" not in ref and not (intro / Path(name).parent / ref).resolve().is_file()]
    assert missing == []


def test_every_file_the_page_sources_refer_to_exists_in_the_repository():
    """The checked-in page refers to modules and files by relative path; each is a file in the repository,
    apart from the generated _assets folder."""
    page = REPO / "web/spike/page"
    missing = [f"{name}: {ref}" for name, ref in cs.referenced_files(str(page))
               if "/_assets/" not in ref and not (page / ref).resolve().is_file()]
    assert missing == []
