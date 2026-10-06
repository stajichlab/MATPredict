"""Aligner comparison harness (read-only; imports the 64971b6 build code exported
to code/). Swaps only the alignment step of classifier_build.build_hmm."""
import hashlib, os, subprocess, sys, tempfile
from pathlib import Path
R = Path(__file__).resolve().parent
CODE = R / "code"
sys.path.insert(0, str(CODE / "src"))
import pyhmmer  # noqa
import MATPredict.detect.classifier_build as cb  # noqa
from MATPredict.detect.family_registry import load_all_families  # noqa

MAFFT = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/basidio-anchors/.pixi/envs/test/bin/mafft"
FAMSA = "/opt/linux/rocky/8.x/x86_64/pkgs/famsa/2.4.1/bin/famsa"
MUSCLE = "/opt/linux/rocky/8.x/x86_64/pkgs/muscle/5.3/bin/muscle"
os.environ["PATH"] = str(Path(MAFFT).parent) + ":" + os.environ["PATH"]

ARMS = {
    "mafft_auto": [MAFFT, "--auto", "--quiet", "--thread", "1"],
    "mafft_linsi": [MAFFT, "--localpair", "--maxiterate", "1000", "--quiet", "--thread", "1"],
    "mafft_einsi": [MAFFT, "--genafpair", "--maxiterate", "1000", "--quiet", "--thread", "1"],
    "muscle5": [MUSCLE],        # -align IN -output OUT
    "famsa": [FAMSA],           # famsa -t 1 IN OUT
}


def align(arm, fa: Path, aln: Path):
    cmd = ARMS[arm]
    if arm.startswith("mafft"):
        with open(aln, "w") as out:
            subprocess.run([*cmd, str(fa)], stdout=out, stderr=subprocess.DEVNULL, check=True)
    elif arm == "muscle5":
        subprocess.run([MUSCLE, "-align", str(fa), "-output", str(aln), "-threads", "1"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    elif arm == "famsa":
        subprocess.run([FAMSA, "-t", "1", str(fa), str(aln)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return aln


def patch(arm, trim=None):
    """Replace cb._mafft with `arm`; optionally trim the alignment (`trim` = callable(aln)->aln)."""
    def _aligner(seqs, workdir, tag):
        fa, aln = workdir / f"{tag}.faa", workdir / f"{tag}.afa"
        cb.write_fasta(fa, cb._sorted(seqs))
        align(arm, fa, aln)
        if trim is not None:
            aln = trim(aln)
        return aln
    cb._mafft = _aligner


def family():
    return next(f for f in load_all_families(CODE / "db")
                if f.key.phylum == "Mucoromycota" and f.key.locus_name == "MAT")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
