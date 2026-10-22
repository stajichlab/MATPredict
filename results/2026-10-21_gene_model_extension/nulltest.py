"""Is the appended 3' region (ext3) homologous to the reference C-terminus? compare its aligned identical columns with shuffled copies of itself."""
import sys, random; sys.path.insert(0, '.')
import modelkit2 as mk, yaml
from pathlib import Path
refs = mk.read_fa("refs_MAT1-1-1.faa")
random.seed(1)
for line in open("glist2.txt"):
    grp, g = line.split()
    if grp not in ("ptr", "teres", "aureo", "macrophomina", "ctrl_alternaria"): continue
    d = yaml.load(open(f"runs/runs190r/{g}/detection_report.yaml"), Loader=mk.L)
    for l in d["detected"]:
        for e in l["gene_evidence"]:
            if e["gene"] != "MAT1-1-1" or not e.get("exons"): continue
            sq = mk.contigs(Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes") / f"{g}.fa.gz", {e["contig"]})[e["contig"]]
            m = mk.gene_models(e, sq, 1)["exo"]
            if m["ext_state"] != "stop_found" or m["ext_aa"] < 30: continue
            ext = m["ext_seq"]
            # pick best of a few refs for the real extension; compare the same ref against shuffles
            best = None
            for rn, rs in refs.items():
                cov, idn, over, nb, ni = mk.align_score(ext, rs)
                if best is None or ni > best[0]: best = (ni, rn, nb)
            rn = best[1]; real = best[0]
            sh = []
            for _ in range(10):
                x = list(ext); random.shuffle(x)
                sh.append(mk.align_score("".join(x), refs[rn])[4])
            print(grp, g[:22], "ext_aa", len(ext), "best ref", rn.split("|")[1][:25], "identical cols real", real, "shuffled mean/max", round(sum(sh)/10,1), max(sh), flush=True)
