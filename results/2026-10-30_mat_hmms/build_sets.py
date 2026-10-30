"""Training sets for a1 / alpha2 / a2 / alpha1 profiles: curated db proteins plus titled NCBI public proteins, one per species per class."""
import re, collections
TAXON = {"4952": ("Yarrowia lipolytica", "Dipodascales"), "1005962": ("Ogataea polymorpha", "Pichiaceae"), "460523": ("Ogataea parapolymorpha", "Pichiaceae"),
         "460519": ("Komagataella phaffii", "Pichiaceae"), "4909": ("Pichia kudriavzevii", "Pichiaceae"), "28985": ("Kluyveromyces lactis", "Saccharomycetaceae"),
         "4932": ("Saccharomyces cerevisiae", "Saccharomycetaceae"), "4929": ("Saccharomyces cerevisiae", "Saccharomycetaceae"), "381046": ("Lachancea kluyveri", "Saccharomycetaceae"),
         "4959": ("Debaryomyces hansenii", "Serinales"), "36911": ("Clavispora lusitaniae", "Serinales"), "5476": ("Candida albicans", "Serinales"),
         "498019": ("Candida auris", "Serinales"), "42374": ("Candida auris", "Serinales")}
SACC = "Monosporozyma Saccharomyces Zygosaccharomyces Kluyveromyces Lachancea Nakaseomyces Kazachstania Tetrapisispora Maudiozyma Naumovozyma Torulaspora Eremothecium Vanderwaltozyma Eremothecium Ashbya Zygotorulaspora Kluyveromyces".split()
SER = "Candidozyma Meyerozyma Candida Clavispora Debaryomyces Meyerozyma Scheffersomyces Lodderomyces Spathaspora Millerozyma Priceomyces Wickerhamiella".split()
PICH = "Hyphopichia Ogataea Komagataella Pichia Kuraishia".split()
def clade(genus):
    if genus in SACC: return "Saccharomycetaceae"
    if genus in PICH: return "Pichiaceae"
    if genus in ("Yarrowia", "Geotrichum", "Dipodascus"): return "Dipodascales"
    if genus in SER: return "Serinales"
    return "other"
import glob, yaml
def species_of(rec):
    f = glob.glob(f"../../db/Ascomycota/*/{rec}/metadata.yaml")[0]
    d = yaml.safe_load(open(f)); sp = d["organism"]["species"]
    return sp, clade(sp.split()[0])
rows = {}   # (class, species) -> (id, seq, clade, source)
def add(cls, species, ident, seq, cl, src):
    k = (cls, species)
    if k not in rows or (src == "curated" and rows[k][3] != "curated") or (src == rows[k][3] and len(seq) > len(rows[k][1])):
        rows[k] = (ident, seq, cl, src)
for l in open("data/curated_mat_proteins.faa").read().split(">")[1:]:
    h, *s = l.split("\n"); rec, name, c = h.split("|"); sp, cl = species_of(rec)
    add(c, sp, f"{rec}|{name}", "".join(s), cl, "curated")
PAT = {"a1": re.compile(r"(mat\s?-?a1|mating[- ]type (protein )?a1|\ba1\b.*transcription|a1p|\bA1 transcription)", re.I),
       "alpha2": re.compile(r"(alpha2|mat\s?-?alpha2)", re.I), "a2": re.compile(r"(mat\s?-?a2|mating[- ]type (protein )?a2|mata2)", re.I)}
for c in ("a1", "alpha2", "a2"):
    for l in open(f"data/public_{c}.faa").read().split(">")[1:]:
        h, *s = l.split("\n"); seq = "".join(s)
        m = re.match(r"(\S+) (.*) \[(.*?)\]\s*$", h)
        if not m or h.startswith("pdb") or "hypothetical" in h or "homeobox" in h.lower() or "partial" in h or "homeodomain" in h.lower() and c == "a2": continue
        acc, title, sp = m.groups(); sp = re.sub(r"\s+(CBS|IFO|NBRC|CD|S288C|DMKU).*$", "", sp.strip("[] "))
        sp = " ".join(sp.replace("[", "").replace("]", "").split()[:2])
        if not PAT[c].search(title) or (c != "alpha2" and re.search(r"alpha", title, re.I)): continue
        if c == "a1" and re.search(r"alpha", title, re.I): continue
        add(c, sp, acc + "|" + title[:40].replace(" ", "_"), seq, clade(sp.split()[0].strip("[]")), "ncbi")
out = collections.defaultdict(list)
with open("data/set_meta.tsv", "w") as m:
    m.write("class\tspecies\tclade\tsource\tid\tlen\n")
    for (c, sp), (i, s, cl, src) in sorted(rows.items()):
        out[c].append((i, s)); m.write(f"{c}\t{sp}\t{cl}\t{src}\t{i}\t{len(s)}\n")
for c, v in out.items():
    open(f"data/set_{c}.faa", "w").write("".join(f">{i.replace(' ', '_')}\n{s}\n" for i, s in v))
    print(c, len(v), collections.Counter(rows[(c, sp)][2] for (cc, sp) in rows if cc == c))
