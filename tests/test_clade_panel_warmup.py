"""The panel script's taxonomy warm-up counts an empty answer as a failure and can stop the job (#324)."""
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_clade_panel.slurm"


def test_the_script_is_valid_bash():
    assert subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True).returncode == 0


def test_an_empty_lineage_and_phylum_counts_as_a_failure():
    text = SCRIPT.read_text()
    assert "if not lineage and not phylum:" in text
    assert "returned no lineage and no phylum" in text
    assert "sys.exit(3 if bad else 0)" in text


def test_strict_mode_is_opt_in_and_defaults_to_a_warning():
    text = SCRIPT.read_text()
    assert 'TAXONOMY_WARMUP_STRICT:-0' in text      # default off: the job still runs, with a warning
    assert '[[ "${TAXONOMY_WARMUP_STRICT:-0}" == "1" ]] && exit 1' in text
