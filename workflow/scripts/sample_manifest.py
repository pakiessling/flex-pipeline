"""Build an auditable manifest before integration, including excluded samples."""
import argparse
import csv
import json
from pathlib import Path


def build_manifest(h5ads, statuses, output):
    if len(h5ads) != len(statuses):
        raise ValueError("Each QC H5AD must have exactly one status file")
    rows = []
    seen = set()
    for path, status_path in zip(h5ads, statuses):
        with open(status_path) as fh:
            status = json.load(fh)
        if status["sample"] in seen:
            raise ValueError(f"Duplicate sample: {status['sample']}")
        seen.add(status["sample"])
        rows.append({
            "sample": status["sample"],
            "h5ad_path": str(Path(path).resolve()),
            "status": "included" if status["integration_eligible"] else "excluded",
            "reason": status["reason"],
            "n_cells_input": status["n_cells_input"],
            "n_cells_retained": status["n_cells_output"],
            "n_cells_integrated": status["n_cells_output"] if status["integration_eligible"] else 0,
            "n_zero_count_cells_removed": status["n_zero_count_cells_removed"],
            "n_variable_genes": status["n_variable_genes"],
            "n_variable_genes_normalized": status["n_variable_genes_normalized"],
        })
    if not rows:
        raise ValueError("No samples were requested")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ads", nargs="+", required=True)
    parser.add_argument("--statuses", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    build_manifest(args.h5ads, args.statuses, args.output)
