#!/usr/bin/env python3
"""Build manifest.tsv / selected.tsv / summary stats for public Fusarium Illumina PE WGS runs.

Data source: ENA portal API (metadata only; no read data is downloaded).
  python3 build_manifest.py CACHE_DIR [--refresh]

CACHE_DIR holds ena_runs.tsv (one ENA search) and tax/<taxid>.json (ENA taxonomy
lookups, used only to assign the species complex from the NCBI lineage).  Missing
files are fetched politely (sleep between calls); --refresh re-fetches everything.
Outputs are written next to this script.
"""
import csv, json, os, re, sys, time, urllib.parse, urllib.request, collections as C

HERE = os.path.dirname(os.path.abspath(__file__))
GENOME_MB = 40.0
LOW_COV = 20.0
FIELDS = ("run_accession,experiment_accession,sample_accession,secondary_sample_accession,"
          "study_accession,scientific_name,tax_id,strain,isolate,library_strategy,library_source,"
          "library_layout,instrument_platform,instrument_model,read_count,base_count,fastq_ftp,"
          "fastq_md5,fastq_bytes,center_name,first_public,country,collection_date,host,"
          "sample_title,study_title,sample_alias")
QUERY = ('tax_tree(5506) AND library_strategy="WGS" AND library_source="GENOMIC" '
         'AND instrument_platform="ILLUMINA" AND library_layout="PAIRED"')

# ---- selection parameters (see summary.md) ---------------------------------
CAP_DEFAULT = 50
CAP_SPECIES = {"Fusarium graminearum": 150, "Fusarium oxysporum": 150}
FSP_CAP = 6            # max runs per named forma specialis inside the F. oxysporum cap
FSP_NONE_CAP = 30      # runs with no forma specialis recorded
COV_SCORE_CEILING = 60.0
# Fola cohort already typed by N. Lee (see summary.md): BioProject PRJNA1098703
USED_EXACT = {"SRR28734943": "yes: 20a (docs/superpowers/specs/2026-10-05-reads-typing...)",
              "SRR28734937": "yes: JCP043 (results/2026-10-05_fola_detect)"}
USED_PROJECT = {"PRJNA1098703": "likely: Fola cohort (SRA runs under nicolel/Fola/02_fastq/NCBI_SRAs_retrieved; accessions not listed in repo)"}

FGSC = {"graminearum", "asiaticum", "boothii", "meridionale", "gerlachii", "aethiopicum",
        "louisianense", "cortaderiae", "vorax", "ussurianum", "mesoamericanum",
        "austroamericanum", "brasilicum", "nepalense", "acaciae-mearnsii"}
LINEAGE_GROUP = {"Fusarium oxysporum species complex": "oxysporum",
                 "Fusarium fujikuroi species complex": "fujikuroi_complex",
                 "Fusarium solani species complex": "solani"}
ENGINEERED_RE = re.compile(r"mutant|crispr|cas9|edit|suppressor|knock|deletion|disrupt|transform|"
                           r"phenotype or genotype|off-target|overexpress|complement|single-asc|"
                           r"TRI5-GFP|\bgfp\b", re.I)


def get(url, dest, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r, open(dest, "wb") as o:
                o.write(r.read())
            return
        except Exception as e:  # noqa
            sys.stderr.write(f"retry {url}: {e}\n"); time.sleep(5)
    raise SystemExit(f"failed: {url}")


def fetch(cache, refresh):
    os.makedirs(os.path.join(cache, "tax"), exist_ok=True)
    runs = os.path.join(cache, "ena_runs.tsv")
    if refresh or not os.path.exists(runs):
        q = urllib.parse.urlencode({"result": "read_run", "query": QUERY, "fields": FIELDS,
                                    "format": "tsv", "limit": 0})
        get("https://www.ebi.ac.uk/ena/portal/api/search?" + q, runs)
    rows = list(csv.DictReader(open(runs), delimiter="\t"))
    for t in sorted({r["tax_id"] for r in rows}):
        f = os.path.join(cache, "tax", t + ".json")
        if refresh or not os.path.exists(f) or os.path.getsize(f) == 0:
            get("https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/" + t.split(";")[0], f)
            time.sleep(0.25)
    return rows


def lineage_group(cache, taxid):
    f = os.path.join(cache, "tax", taxid + ".json")
    try:
        j = json.load(open(f)); j = j[0] if isinstance(j, list) else j
        for part in j.get("lineage", "").split("; "):
            if part in LINEAGE_GROUP:
                return LINEAGE_GROUP[part], part
        return "", j.get("lineage", "")
    except Exception:
        return "", "taxonomy_lookup_failed"


def species_of(name):
    n = name.split(";")[0].strip()
    m = re.match(r"(Fusarium\s+(?:cf\.\s+|aff\.\s+)?[a-z\-]+)", n)
    if not m: return n
    return "Fusarium sp." if m.group(1) == "Fusarium sp" else m.group(1)


def fsp_of(r):
    for field in ("scientific_name", "strain", "isolate", "sample_title", "sample_alias", "study_title"):
        m = re.search(r"f\.\s*sp\.\s*([A-Za-z\-]+)", r[field])
        if m:
            return m.group(1).lower(), field
    return "", ""


def ctry(o):
    """Country used for spreading: text before ':'; missing-value strings become 'unknown'."""
    c = o["country"].split(":")[0].strip()
    return "unknown" if (not c or re.match(r"(?i)missing|not |na$|n/a|unknown|none", c)) else c


def pair_files(r):
    """Return (urls, md5s, bytes, extra_note) restricted to _1/_2 mates."""
    urls = r["fastq_ftp"].split(";") if r["fastq_ftp"] else []
    md5 = r["fastq_md5"].split(";") if r["fastq_md5"] else []
    sz = r["fastq_bytes"].split(";") if r["fastq_bytes"] else []
    if len(md5) != len(urls) or len(sz) != len(urls):
        return [], [], [], "ftp/md5/bytes_count_mismatch"
    idx = {}
    for i, u in enumerate(urls):
        m = re.search(r"_([12])\.fastq\.gz$", u)
        if m:
            idx[m.group(1)] = i
    if set(idx) != {"1", "2"}:
        return [], [], [], f"no_mate_pair_in_ENA_files({len(urls)}_files)"
    sel = [idx["1"], idx["2"]]
    note = "extra_unpaired_file_ignored" if len(urls) > 2 else ""
    return [urls[i] for i in sel], [md5[i] for i in sel], [sz[i] for i in sel], note


def main():
    cache = sys.argv[1]
    rows = fetch(cache, "--refresh" in sys.argv)
    epithets = set()
    for r in rows:
        m = re.match(r"Fusarium\s+(?:cf\.\s+|aff\.\s+)?([a-z\-]+)", r["scientific_name"])
        if m: epithets.add(m.group(1))
    epithets -= {"sp", "oxysporum"}  # 'oxysporum' handled explicitly below

    out = []
    for r in rows:
        sp = species_of(r["scientific_name"])
        epi = sp.split()[-1]
        pooled = ";" in r["tax_id"]
        if pooled: sp = "Fusarium oxysporum"  # all pooled names are F. oxysporum strains/unidentified
        grp, cplx = lineage_group(cache, r["tax_id"])
        if pooled:
            grp = grp or "oxysporum"
        if not grp:
            grp = "graminearum_complex" if epi in FGSC and sp.startswith("Fusarium ") else "other"
        elif epi in FGSC:
            grp = "graminearum_complex"
        fsp, fsp_src = fsp_of(r)
        urls, md5, sz, pnote = pair_files(r)
        bc = int(r["base_count"]) if r["base_count"] else 0
        notes = []
        if pnote: notes.append(pnote)
        if pooled: notes.append("pooled_multi_isolate_run(tax_id list)")
        if fsp and fsp_src != "scientific_name": notes.append(f"f.sp._from_{fsp_src}")
        text = " ".join(r[k] for k in ("strain", "isolate", "sample_title", "sample_alias"))  # sample-level only
        if ENGINEERED_RE.search(r["study_title"] + " " + r["sample_title"] + " " + r["strain"] + " " + r["isolate"]):
            notes.append("engineered_or_lab_mutant_suspect(keyword in study/sample title)")
        # label suspect
        sus = []
        if re.search(r"\b(cf\.|aff\.)", r["scientific_name"]) or sp in ("Fusarium sp.", "unidentified"):
            sus.append("uncertain_label(cf./sp./unidentified)")
        for m in re.finditer(r"\b(?:Fusarium|F\.)\s+([a-z]{4,})", text):
            e = m.group(1)
            if e != epi and (e in epithets or e == "oxysporum") and e not in ("species", "genome"):
                sus.append(f"text_mentions_{e}")
        sus = sorted(set(sus))
        used = USED_EXACT.get(r["run_accession"]) or USED_PROJECT.get(r["study_accession"], "")
        out.append(dict(
            run_accession=r["run_accession"], sample_accession=r["secondary_sample_accession"],
            biosample=r["sample_accession"], bioproject=r["study_accession"], species=sp,
            species_group=grp, taxid=r["tax_id"] if not pooled else "multiple", scientific_name_ena=r["scientific_name"][:80],
            strain_isolate=(r["strain"] or r["isolate"] or r["sample_alias"]),
            forma_specialis=fsp, country=r["country"], host=r["host"], center=r["center_name"],
            instrument_model=r["instrument_model"], read_count=r["read_count"], base_count=r["base_count"],
            est_coverage=f"{bc / (GENOME_MB * 1e6):.1f}", fastq_ftp=";".join(urls), fastq_bytes=";".join(sz),
            md5=";".join(md5), first_public=r["first_public"], collection_date=r["collection_date"],
            low_coverage="Y" if bc / (GENOME_MB * 1e6) < LOW_COV else "",
            duplicate_sample="", species_label_suspect=";".join(sus),
            already_used_in_matpredict=used, notes=";".join(notes),
            _cov=bc / (GENOME_MB * 1e6), _pairs=bool(urls), _pooled=pooled,
            _eng=any(n.startswith("engineered") for n in notes), _study_title=r["study_title"]))
    # duplicates: several runs of one biosample; and same strain under several biosamples
    bs = C.Counter(o["biosample"] for o in out)
    st = C.defaultdict(set)
    for o in out:
        if o["strain_isolate"] and not o["_pooled"]:
            st[(o["species"], o["strain_isolate"].strip().lower())].add(o["biosample"])
    for o in out:
        if bs[o["biosample"]] > 1: o["duplicate_sample"] = f"Y({bs[o['biosample']]}_runs_in_biosample)"
        k = (o["species"], o["strain_isolate"].strip().lower())
        if o["strain_isolate"] and len(st.get(k, ())) > 1:
            o["notes"] = (o["notes"] + ";" if o["notes"] else "") + f"same_strain_name_in_{len(st[k])}_biosamples"
    out.sort(key=lambda o: (o["species_group"], o["species"], o["bioproject"], o["run_accession"]))
    cols = [c for c in out[0] if not c.startswith("_")]
    with open(os.path.join(HERE, "manifest.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(cols)
        for o in out: w.writerow([o[c] for c in cols])

    # ---- selection ----------------------------------------------------------
    elig = [o for o in out if o["_pairs"] and not o["_pooled"] and not o["_eng"]
            and o["low_coverage"] == "" and o["already_used_in_matpredict"] == ""
            and o["bioproject"] != "PRJNA1098703"]
    # one run per biosample: highest coverage
    best = {}
    for o in elig:
        b = best.get(o["biosample"])
        if b is None or o["_cov"] > b["_cov"]: best[o["biosample"]] = o
    cand = C.defaultdict(list)
    for o in best.values(): cand[o["species"]].append(o)
    sel = []
    for spc, lst in cand.items():
        cap = CAP_SPECIES.get(spc, CAP_DEFAULT)
        chosen = []; nbp = C.Counter(); nct = C.Counter(); nfs = C.Counter()
        pool = list(lst)
        while pool and len(chosen) < cap:
            def key(o):
                fs = o["forma_specialis"] or "(none)"
                return ((nfs[fs] if spc == "Fusarium oxysporum" else 0),
                        nbp[o["bioproject"]], nct[ctry(o)],
                        -min(o["_cov"], COV_SCORE_CEILING), sum(int(x) for x in o["fastq_bytes"].split(";")))
            if spc == "Fusarium oxysporum":
                pool = [o for o in pool if nfs[o["forma_specialis"] or "(none)"] < (FSP_CAP if o["forma_specialis"] else FSP_NONE_CAP)]
                if not pool: break
            o = min(pool, key=key); pool.remove(o); chosen.append(o)
            o["selection_rank"] = len(chosen)  # order picked within species (truncate to shrink)
            nbp[o["bioproject"]] += 1; nct[ctry(o)] += 1; nfs[o["forma_specialis"] or "(none)"] += 1
        sel += chosen
    sel.sort(key=lambda o: (o["species_group"], o["species"], o["run_accession"]))
    scols = cols + ["selection_rank"]
    with open(os.path.join(HERE, "selected.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(scols)
        for o in sel: w.writerow([o[c] for c in scols])
    print(len(out), "manifest rows;", len(elig), "eligible runs;", len(best), "eligible biosamples;", len(sel), "selected")


if __name__ == "__main__":
    main()
