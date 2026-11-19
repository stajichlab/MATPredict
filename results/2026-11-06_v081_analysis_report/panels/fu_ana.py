import yaml,glob,os,collections,statistics as st
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"; F=R+"/2026-11-06_seg_followups"
def rd(p):
    try: return yaml.safe_load(open(p+"/detection_report.yaml"))
    except Exception: return None
def loci(r): return sorted((x["family"],x["idiomorph"],x["confidence"],x["contig"],x["start"],x["end"]) for x in (r.get("detected") or []))
def coarse(r): return sorted((x["family"],x["idiomorph"],x["contig"]) for x in (r.get("detected") or []))
def genes(r): return sorted((x["family"],x["contig"],tuple(sorted({e["gene"] for e in x["gene_evidence"] if not e.get("lost_to")}))) for x in (r.get("detected") or []))
def wall(p):
    try: return int(open(p+"/wall_seconds").read())
    except Exception: return None
def compare(name,pairs,strata):
    n=same=cs=lost=gm=0; tr=[]; byst=collections.defaultdict(lambda:[0,0]); diffs=[]; missing=0
    for g,(bp,np_) in pairs.items():
        a,b=rd(bp),rd(np_)
        if not a or not b: missing+=1; continue
        n+=1; s=strata.get(g,"?"); byst[s][1]+=1
        if loci(a)==loci(b): same+=1; byst[s][0]+=1
        else:
            if coarse(a)==coarse(b): cs+=1
            else: lost+=1; diffs.append((g,s,coarse(a),coarse(b)))
        if genes(a)!=genes(b): gm+=1
        wa,wb=wall(bp),wall(np_)
        if wa and wb: tr.append((wb,wa))
    print(f"\n=== {name}: compared {n}, missing/failed {missing}")
    print(f"identical loci {same} ({100*same/max(n,1):.1f}%); same family/idiomorph/contig but coordinates/confidence differ {cs}; family/idiomorph/contig differ {lost}; gene set of a locus differs in {gm} genomes")
    print("per stratum identical/total:",{k:tuple(v) for k,v in sorted(byst.items())})
    if tr: print("wall seg-on / baseline: median %.2f, total %d s vs %d s (%.2fx)"%(st.median(a/b for a,b in tr),sum(a for a,_ in tr),sum(b for _,b in tr),sum(a for a,_ in tr)/sum(b for _,b in tr)))
    for g,s,a,b in diffs[:12]: print("  DIFF",g[:40],s,"baseline",a,"| new",b)
# Basidiomycota
strata={l.split("\t")[0]:l.rstrip().split("\t")[2] for l in open(F+"/bas_panel.tsv")}
pairs={}
for g in strata:
    bp=(glob.glob(f"{R}/2026-11-05_v081_basidiomycota/runs/*/runs/{g}")+[None])[0]; npth=f"{F}/bas/runs/{g}"
    if bp and os.path.isdir(npth): pairs[g]=(bp,npth)
compare("Basidiomycota (SEG on + soft masking vs v0.8.1)",pairs,strata)
# Mucoromycotina
strata={l.split("\t")[0]:l.rstrip().split("\t")[2] for l in open(F+"/muc_panel.tsv")}
pairs={g:(f"{R}/2026-11-05_v081_mucoromycotina/runs/{g}",f"{F}/muc/{g}") for g in strata if os.path.isdir(f"{F}/muc/{g}")}
compare("Mucoromycotina (SEG on + soft masking vs v0.8.1)",pairs,strata)
# same-node timing, 60 Ascomycota genomes
base={l.split()[0]:l.split()[1] for l in open(R+"/2026-11-06_seg_call_test/baseline_paths.txt")}
ids=[l.split("\t")[0] for l in open(F+"/lists/asco60.tsv")]
tn={g:f"{F}/t_no/runs/{g}" for g in ids if os.path.isdir(f"{F}/t_no/runs/{g}")}; ty={g:f"{F}/t_yes/runs/{g}" for g in ids if os.path.isdir(f"{F}/t_yes/runs/{g}")}
common=[g for g in ids if g in tn and g in ty and wall(tn[g]) and wall(ty[g])]
print(f"\n=== same-node timing (r12, 16 genomes at a time), {len(common)} genomes")
if common:
    a=[wall(tn[g]) for g in common]; b=[wall(ty[g]) for g in common]
    print("SEG off (production): median %.0f s, total %d s | SEG on + soft masking: median %.0f s, total %d s | ratio of totals %.2f, median per-genome ratio %.2f"%(st.median(a),sum(a),st.median(b),sum(b),sum(b)/sum(a),st.median(y/x for x,y in zip(a,b))))
    slow=sorted(common,key=lambda g:-wall(tn[g]))[:5]; print("slowest five (off -> on):",[(g[:22],wall(tn[g]),wall(ty[g])) for g in slow])
compare("same-node SEG off run vs v0.8.1 baseline (determinism)",{g:(base[g],tn[g]) for g in tn if g in base},{})
compare("same-node SEG on run vs SEG off run",{g:(tn[g],ty[g]) for g in common},{})
