#!/usr/bin/env bash
# launch.sh — Start the flex-pipeline on a SLURM HPC cluster.
#
# Usage:
#   ./launch.sh                    # run with defaults from config/config.yaml
#   ./launch.sh --dry-run          # preview jobs without submitting
#   ./launch.sh --local            # run locally (no SLURM, useful for testing)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$SCRIPT_DIR/config/config.yaml"
CLUSTER_CFG="$SCRIPT_DIR/config/cluster.yaml"

# ---------------------------------------------------------------------------
# Validate prerequisites
# ---------------------------------------------------------------------------
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: $CONFIG not found. Have you set up the pipeline?" >&2
    exit 1
fi

SAMPLES_CSV=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG')).get('samples','samples.csv'))")
if [ ! -f "$SAMPLES_CSV" ]; then
    echo "ERROR: samples CSV not found: $SAMPLES_CSV" >&2
    echo "Create it with columns: sample_id, raw_matrix_path, filtered_matrix_path" >&2
    exit 1
fi

if ! command -v snakemake &>/dev/null; then
    echo "ERROR: snakemake not found in PATH." >&2
    echo "Activate the conda environment that has Snakemake installed, or set snakemake_env in config/cluster.yaml." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Parse cluster config
# ---------------------------------------------------------------------------
_yaml_val() {
    python3 -c "
import yaml, sys
d = yaml.safe_load(open('$CLUSTER_CFG')) if __import__('os').path.exists('$CLUSTER_CFG') else {}
print(d.get('$1', ''))
" 2>/dev/null || true
}

ACCOUNT=$(_yaml_val slurm_account)
PARTITION=$(_yaml_val slurm_partition)

# ---------------------------------------------------------------------------
# Parse flags
# ---------------------------------------------------------------------------
DRY_RUN=""
LOCAL=""
for arg in "$@"; do
    case "$arg" in
        --dry-run)  DRY_RUN="--dry-run" ;;
        --local)    LOCAL="1" ;;
        *)          echo "Unknown option: $arg" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Build snakemake command
# ---------------------------------------------------------------------------
if [ -n "$LOCAL" ]; then
    CMD="snakemake \
        --snakefile workflow/Snakefile \
        --configfile config/config.yaml \
        --conda-frontend conda \
        --use-conda \
        --cores 4 \
        $DRY_RUN"
else
    # Generate a runtime profile that merges profiles/config.yaml with the
    # cluster-specific account and partition from config/cluster.yaml.
    # This avoids the snakemake profile vs CLI --default-resources conflict,
    # where the profile's default-resources block silently wins over CLI flags.
    RUNTIME_PROFILE_DIR="$SCRIPT_DIR/.runtime_profile"
    mkdir -p "$RUNTIME_PROFILE_DIR"
    python3 -c "
import yaml
with open('$SCRIPT_DIR/profiles/config.yaml') as f:
    profile = yaml.safe_load(f)
dr = profile.setdefault('default-resources', {})
dr['slurm_account'] = '$ACCOUNT'
dr['slurm_partition'] = '$PARTITION'
with open('$RUNTIME_PROFILE_DIR/config.yaml', 'w') as f:
    yaml.dump(profile, f, default_flow_style=False)
"
    echo "Runtime profile written to $RUNTIME_PROFILE_DIR/config.yaml"
    echo "  slurm_account:   $ACCOUNT"
    echo "  slurm_partition: $PARTITION"
    echo ""

    CMD="snakemake \
        --snakefile workflow/Snakefile \
        --configfile config/config.yaml \
        --executor slurm \
        --conda-frontend conda \
        --profile $RUNTIME_PROFILE_DIR \
        --jobscript jobscript.sh \
        --rerun-incomplete \
        $DRY_RUN"
fi

echo "Running: $CMD"
echo ""

if [ -n "$DRY_RUN" ]; then
    eval "$CMD"
else
    mkdir -p logs
    nohup bash -c "$CMD" > logs/pipeline_$(date +%Y%m%d_%H%M%S).log 2>&1 &
    echo "Pipeline started in background. PID: $!"
    echo "Monitor progress: tail -f logs/pipeline_*.log"
fi
