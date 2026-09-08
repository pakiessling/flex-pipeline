from pathlib import Path
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("is_log1p, flag", [(False, "--no-is_log1p"), (True, "--is_log1p")])
def test_markers_only_target_and_boolean_flag(tmp_path, is_log1p, flag):
    source = tmp_path / "input.h5ad"
    source.touch()
    samples = tmp_path / "samples.csv"
    samples.write_text(f"sample_id,h5ad_path\nA,{source}\n")
    config = {
        "samples": str(samples),
        "steps": {"qc": True, "qc_diagnostics": False, "integration": True,
                  "integration_diagnostics": False, "markers": True,
                  "label_transfer": False, "cytetype": False, "report": False},
        "params": {"leiden_resolutions": [1.0], "leiden_primary_resolution": 1.0,
                   "markers_is_log1p": is_log1p},
        "output": {key: str(tmp_path / key) for key in (
            "intermediate", "per_sample", "qc_plots", "integration",
            "annotation", "markers", "report")},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    result = subprocess.run(
        [sys.executable, "-m", "snakemake", "--snakefile", str(ROOT / "workflow/Snakefile"),
         "--configfile", str(config_path), "--cores", "1", "--dry-run", "--printshellcmds",
         "--shared-fs-usage", "persistence", "software-deployment", "sources", "input-output"],
        cwd=ROOT, capture_output=True, text=True,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "marker_genes_leiden_1_0.csv" in output
    assert flag in output
    if is_log1p:
        assert "--no-is_log1p" not in output
