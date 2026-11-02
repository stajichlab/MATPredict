#!/usr/bin/env python3
"""Check an assembled site (the output of build-site.sh) before it is uploaded or published.

It catches what a green build does not: a page module that was never copied into the site, a file the page
fetches that is not there, and a database archive that does not match its manifest. Each of those has
shipped, or nearly shipped, a page that loads and then fails when a genome is run.

Usage: check_site.py <site-dir>        (exit 0 when clean, 1 with one line per problem otherwise)
"""
import hashlib
import json
import os
import re
import sys

REQUIRED = (
    "index.html", "results.html", "about.html", "data/summary.json",
    "page/index.html", "page/app.mjs",
    "page/_assets/matpredict.whl", "page/_assets/db-manifest.json", "page/_assets/db.tar.gz",
)

# A quoted relative path with a file extension: an import ("./genbank.mjs", "../browser/host/tool-host.mjs") or a
# fetched file ("../localizer/localizer.py"). Template strings (`${ASSETS}/...`) are not matched.
REFERENCE = re.compile(r"""["'](\.{1,2}/[^"'\s${}]+\.[A-Za-z0-9]+)["']""")


def referenced_files(page_dir: str) -> list[tuple[str, str]]:
    """(module name, relative path) for every quoted relative file reference in the page's .mjs files."""
    out = []
    for name in sorted(os.listdir(page_dir)):
        if name.endswith(".mjs"):
            with open(os.path.join(page_dir, name), encoding="utf-8") as fh:
                for ref in REFERENCE.findall(fh.read()):
                    out.append((name, ref))
    return out


def check(site: str) -> list[str]:
    problems = []
    for rel in REQUIRED:
        if not os.path.isfile(os.path.join(site, rel)):
            problems.append(f"missing {rel}")
    page = os.path.join(site, "page")
    if os.path.isdir(page):
        for name, ref in referenced_files(page):
            if not os.path.isfile(os.path.normpath(os.path.join(page, ref))):
                problems.append(f"page/{name} refers to {ref}, which is not in the site")
    assets = os.path.join(page, "_assets")
    manifest_path = os.path.join(assets, "db-manifest.json")
    if os.path.isfile(manifest_path):
        manifest = json.load(open(manifest_path))
        archive = manifest.get("dbArchive")
        if archive:
            path = os.path.join(assets, archive)
            if not os.path.isfile(path):
                problems.append(f"the manifest names {archive}, which is not in the site")
            else:
                blob = open(path, "rb").read()
                if len(blob) != manifest.get("dbArchiveBytes"):
                    problems.append(f"{archive} is {len(blob)} bytes, the manifest says {manifest.get('dbArchiveBytes')}")
                if hashlib.sha256(blob).hexdigest() != manifest.get("dbArchiveSha256"):
                    problems.append(f"{archive} does not match the SHA-256 in the manifest")
        wheel = manifest.get("pyhmmerWheel")
        if wheel and not os.path.isfile(os.path.join(assets, "bundle", wheel)):
            problems.append(f"the manifest names the pyhmmer wheel {wheel}, which is not in the bundle")
    return problems


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    problems = check(argv[1])
    for p in problems:
        print(f"site check: {p}")
    if problems:
        return 1
    n = len(referenced_files(os.path.join(argv[1], "page")))
    print(f"site check: ok ({len(REQUIRED)} required files, {n} page references, database archive matches its manifest)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
