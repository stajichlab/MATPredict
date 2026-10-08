"""Pull a clean core gene set (flanks + MAT genes) for one contig from a MATPredict detected_loci.gff3."""
CORE = {"COX13", "APN2", "SLA2", "MAT1-1-1", "MAT1-2-1"}


def gff_genes(gff, contig_prefix, min_identity=50.0):
    best = {}
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or not f[0].startswith(contig_prefix):
            continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        name = a.get("Name")
        ident = float(a.get("identity", 0) or 0)
        if name in CORE and a.get("present") == "true" and ident >= min_identity:
            if name not in best or ident > best[name][4]:
                best[name] = (int(f[3]), int(f[4]), f[6], name, ident)
    return [(s, e, st, n) for (s, e, st, n, _i) in best.values()]


def locus_contig(report_yaml, idiomorph=None):
    """Contig (full name) and interval of the longest detected locus in a detection_report.yaml."""
    import yaml
    y = yaml.safe_load(open(report_yaml))
    loci = [x for x in (y.get("detected") or []) if idiomorph is None or x["idiomorph"] == idiomorph]
    if not loci:
        return None
    x = max(loci, key=lambda z: z["end"] - z["start"])
    return x["contig"], x["start"], x["end"], x["idiomorph"]
