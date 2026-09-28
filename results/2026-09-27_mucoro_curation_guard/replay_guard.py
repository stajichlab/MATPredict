"""Replay the secondary-undetermined guard on existing reports.
Guard: in an enum-vocabulary family, a call with idiomorph 'undetermined' (not
homothallic_candidate) is withheld when the same family has a determined call
in the same genome."""
import glob, io, sys, tarfile, collections, yaml
import subprocess
DB = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/curation-umbelopsis/db"
enum = set()
for f in glob.glob(f"{DB}/*/order.yml"):
    phy = f.split("/")[-2]
    d = yaml.safe_load(open(f))
    for loc in d.get("loci", []):
        if loc.get("vocabulary_type") == "enum":
            enum.add(f"{phy}:{loc['locus_name']}")
def reports(path):
    for f in glob.glob(f"{path}/**/detection_report.yaml", recursive=True):
        yield f, yaml.safe_load(open(f))
    for tz in glob.glob(f"{path}/**/reports*.tar.zst", recursive=True):
        if glob.glob(f"{path}/**/detection_report.yaml", recursive=True):
            break
        data = subprocess.run(["zstd", "-dc", tz], capture_output=True, check=True).stdout
        with tarfile.open(fileobj=io.BytesIO(data)) as t:
            for m in t.getmembers():
                if m.name.endswith("detection_report.yaml"):
                    yield tz + ":" + m.name, yaml.safe_load(t.extractfile(m))
for path in sys.argv[1:]:
    n_g = n_calls = 0; hit = []
    for f, r in reports(path):
        n_g += 1
        det = r.get("detected") or []
        n_calls += len(det)
        by = collections.defaultdict(list)
        for x in det:
            by[x["family"]].append(x)
        for fam, xs in by.items():
            if fam not in enum:
                continue
            determined = [x for x in xs if x.get("idiomorph") not in (None, "undetermined")]
            if not determined:
                continue
            for x in xs:
                if x.get("idiomorph") == "undetermined" and x.get("locus_class") != "homothallic_candidate":
                    hit.append((f.split("/")[-2] if ":" not in f else f.split(":")[1].split("/")[-2], fam, x["contig"], x.get("confidence"), x.get("locus_class"), ",".join(x.get("genes_found", []))))
    print(f"== {path}: genomes {n_g}, calls {n_calls}, withheld by guard {len(hit)}")
    for h in hit[:40]:
        print("   ", *h)
