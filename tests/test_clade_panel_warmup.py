"""The panel script runs offline by default, records the taxonomy snapshot and stops on a taxid with no lineage (#324, #342)."""
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_clade_panel.slurm"


def test_the_script_is_valid_bash():
    assert subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True).returncode == 0


def test_the_warm_up_is_the_tested_helper():
    text = SCRIPT.read_text()
    assert '"$REPO/scripts/taxonomy_warmup.py"' in text     # behaviour tested in tests/test_taxonomy_warmup.py
    assert "if not lineage and not phylum:" not in text      # no second inline copy to drift


def test_offline_is_the_default_and_needs_the_table():
    text = SCRIPT.read_text()
    assert 'MATPREDICT_OFFLINE="${MATPREDICT_OFFLINE:-1}"' in text
    assert "is not a file on" in text and "exit 1" in text      # offline with no table stops the job


def test_strict_warm_up_defaults_on_offline_and_off_online():
    text = SCRIPT.read_text()
    assert 'TAXONOMY_WARMUP_STRICT="${TAXONOMY_WARMUP_STRICT:-1}"' in text
    assert 'TAXONOMY_WARMUP_STRICT="${TAXONOMY_WARMUP_STRICT:-0}"' in text
    assert '[[ "$TAXONOMY_WARMUP_STRICT" == "1" ]] && exit 1' in text


def test_the_snapshot_is_written_to_the_run_folder():
    text = SCRIPT.read_text()
    assert '--describe > "$OUT_DIR/TAXONOMY.txt"' in text


def test_the_offline_block_runs_in_bash_and_stops_without_a_table(tmp_path):
    """Run the offline block alone: with MATPREDICT_OFFLINE on and no table the script exits 1 with a message."""
    text = SCRIPT.read_text()
    start = text.index("export MATPREDICT_OFFLINE=")
    end = text.index("[[ -f \"$REPO/scripts/taxonomy_warmup.py\" ]]")
    block = text[start:end]
    r = subprocess.run(["bash", "-c", block + "\necho passed"], capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin", "MATPREDICT_TAXONOMY": str(tmp_path / "missing.tsv.zst")})
    assert r.returncode == 1 and "is not a file on" in r.stderr and "passed" not in r.stdout
    table = tmp_path / "t.tsv.zst"
    table.write_text("x")
    r = subprocess.run(["bash", "-c", block + '\necho "strict=$TAXONOMY_WARMUP_STRICT"'], capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin", "MATPREDICT_TAXONOMY": str(table)})
    assert r.returncode == 0 and "strict=1" in r.stdout
    r = subprocess.run(["bash", "-c", block + '\necho "strict=$TAXONOMY_WARMUP_STRICT"'], capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin", "MATPREDICT_OFFLINE": "0"})
    assert r.returncode == 0 and "strict=0" in r.stdout
