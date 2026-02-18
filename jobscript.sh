#!/usr/bin/env bash
# jobscript.sh — Snakemake SLURM job script template.
# This file is used by snakemake-executor-plugin-slurm when submitting jobs.
# Do NOT run directly; it is invoked automatically by Snakemake.
set -euo pipefail

# Locate conda root (reads from cluster.yaml if yq is available, else auto-detects)
_detect_conda() {
    local cfg="config/cluster.yaml"
    if command -v python3 &>/dev/null && [ -f "$cfg" ]; then
        local root
        root=$(python3 -c "
import yaml, sys
d = yaml.safe_load(open('$cfg'))
root = d.get('conda_root', '').strip()
print(root)
" 2>/dev/null || true)
        if [ -n "$root" ] && [ -d "$root" ]; then
            echo "$root"
            return
        fi
    fi
    if [ -d "$HOME/miniconda3" ]; then
        echo "$HOME/miniconda3"
    elif [ -d "$HOME/anaconda3" ]; then
        echo "$HOME/anaconda3"
    else
        echo ""
    fi
}

CONDA_ROOT=$(_detect_conda)
if [ -z "$CONDA_ROOT" ]; then
    echo "ERROR: Could not locate conda. Set conda_root in config/cluster.yaml." >&2
    exit 1
fi

# shellcheck disable=SC1091
. "$CONDA_ROOT/etc/profile.d/conda.sh"
export PATH="$CONDA_ROOT/bin:$PATH"

{exec_job}
