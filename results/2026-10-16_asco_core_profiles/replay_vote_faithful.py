"""Replay the idiomorph vote with the pipeline's own supersession, correcting replay_vote.py.

replay_vote.py let every gene vote, but the pipeline first resolves overlapping mutually exclusive pairs
(`idiomorph_resolutions`: the loser is superseded and does not vote). usage: replay_vote_faithful.py REPORTS_DIR
scored.tsv order.yml
Variants: A the campaign rule (losers excluded); B A with optional genes also removed from the vote; C B plus a
resolution step that skips any pair containing an optional gene, so an accessory gene cannot supersede a core hit.
"""
import csv, glob, os, re, sys, collections, yaml
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
sc = {(r["genome"], r["contig"], r["start"]): r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
opt = set(re.findall(r"name:\s*([A-Za-z0-9_.-]+),[^}]*optional:\s*true", open(sys.argv[3]).read()))
def idm(n):
    m = re.match(r"MAT1-([12])-", n); return None if not m else "MAT1-" + m.group(1)
def vote(ev, excluded):
    best = {}
    for e in ev:
        i = idm(e["gene"])
        if not i or e["gene"] in excluded: continue
        s = e.get("bitscore") if e.get("bitscore") is not None else e.get("identity")
        if s is not None and s > best.get(i, -1): best[i] = s
    if not best: return None
    r = sorted(best.items(), key=lambda kv: -kv[1])
    return r[0][0] if len(r) == 1 or r[0][1] > r[1][1] else "undetermined"
tot = collections.Counter(); ch = {v: collections.Counter() for v in "BC"}; prof = {v: collections.Counter() for v in "ABC"}
for p in glob.glob(sys.argv[1] + "/**/detection_report.yaml", recursive=True):
    g = os.path.basename(os.path.dirname(p)); doc = yaml.load(open(p), Loader=L)
    for c in doc["detected"]:
        if c["family"] != "Ascomycota:MAT": continue
        ev = c["gene_evidence"]; res = c.get("idiomorph_resolutions") or []
        losers = {r["loser"] for r in res}
        keep_pair = {r["loser"] for r in res if r["winner"] not in opt and r["loser"] not in opt}
        A = vote(ev, losers); B = vote(ev, losers | opt); Cc = vote(ev, keep_pair | opt)
        tot[A == c["idiomorph"]] += 1
        for name, v in (("B", B), ("C", Cc)):
            if v is not None and v != c["idiomorph"]: ch[name][(c["idiomorph"], v)] += 1
        s = sc.get((g, c["contig"], str(c["start"])))
        if s:
            for name, v in (("A", A), ("B", B), ("C", Cc)):
                if v is not None: prof[name][v == s["profile_call"]] += 1
print("optional genes:", sorted(opt)); print("replay A equals the campaign label:", dict(tot))
for n in "BC": print(f"variant {n}: labels changed {sum(ch[n].values())} {dict(ch[n])}")
for n in "ABC": print(f"variant {n}: agrees with the profile on {prof[n][True]} of {sum(prof[n].values())} scored loci")
