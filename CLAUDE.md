# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a **Snakemake-based bioinformatics pipeline** for processing 10X Genomics FLEX single-cell RNA-seq data. It is organized into two sequential phases:

1. **`01_cleaning/`** — Per-sample QC: ambient RNA removal (SoupX) and quality control filtering
2. **`02_integration/`** — Multi-sample integration: Harmony batch correction, cell type annotation (CellTypist), and marker gene identification

The pipeline runs on SLURM HPC clusters using Conda-managed environments.

## Running the Pipeline

**Phase 1 (per-sample cleaning):**
```bash
cd 01_cleaning
./launch.sh
# Equivalent: snakemake --conda-frontend conda --profile slurm --jobscript jobscript.sh --config samples="samples.csv"
```

**Phase 2 (integration):**
```bash
cd 02_integration
./launch.sh
# Equivalent: snakemake --executor slurm --conda-frontend conda --config input="00_cleaning/output"
```

**Dry-run (preview without execution):**
```bash
snakemake -n --config samples="samples.csv"
```

**Run a single rule locally:**
```bash
snakemake --use-conda <rule_name> --config samples="samples.csv"
```

## Architecture

### Data Flow

```
CellRanger output (raw + filtered H5 matrices)
        ↓ [01_SoupX.py]   - Ambient RNA correction via R/SoupX (rpy2 bridge)
        ↓ [02_qc.py]      - QC metrics, doublet detection, clustering per sample
        ↓ H5AD files (one per sample, in 01_cleaning/output/)
        ↓ [03_integration.py] - Merge samples, Harmony batch correction, UMAP/PaCMAP
        ↓ [04_celltyping.py]  - CellTypist annotation (immune + heart models)
        ↓ [05_markers.R]      - Wilcoxon rank-sum marker genes via Presto
```

### Key Design Decisions

- **AnnData layers** store multiple processing states: `b4_soupx` and `after_soupx`. `03_integration.py` checks which layer to use based on whether SoupX ran successfully.
- **Cells are marked but not removed** during QC (`02_qc.py`). Doublets and low-quality cells are annotated in `.obs` metadata; hard filtering happens in `03_integration.py` (removes `predicted_doublet == True`).
- **SoupX graceful fallback**: if SoupX fails, `01_SoupX.py` falls back to the unmodified filtered matrix rather than crashing the pipeline.
- **Reproducibility**: all scripts set `np.random.seed(0)`, `random.seed(0)`, and `PYTHONHASHSEED=0`.
- **Cell type annotation is sequential**: immune model (`Immune_All_Low.pkl`) runs first, then heart-specific (`Healthy_Adult_Heart.pkl`). Multiple prediction columns with confidence scores are stored in `.obs`.

### Conda Environments

| File | Used by | Key packages |
|------|---------|--------------|
| `01_cleaning/enviroments/py_r.yml` | `run_soupx` rule | Python 3.11, R 4.3.3, SoupX, rpy2, anndata2ri |
| `01_cleaning/enviroments/sc.yml` | `run_qc` rule | ScanPy, Celltypist, DoubletDetection, PaCMAP |
| `02_integration/enviroment.yml` | all phase 2 rules | Harmony, CellTypist, anndataR, Presto, DESeq2 |

### SLURM Configuration

- Phase 1 (`01_cleaning/slurm/config.yaml`): 50 GB RAM, 24hr walltime, max 100 jobs
- Phase 2 (`02_integration/profiles/config.yaml`): 180 GB RAM, 12 CPUs, 12hr walltime, max 500 jobs
- Account: `p0020567`; conda assumed at `$HOME/miniconda3` with fallback to `/work/mz637064/miniconda3/`

### Input Format

Phase 1 requires a `samples.csv` with columns:
- `sample_id` — unique identifier
- `raw_matrix_path` — path to CellRanger raw matrix directory
- `filtered_matrix_path` — path to CellRanger filtered matrix directory

Phase 2 auto-discovers all `.h5ad` files in the configured input directory.

### Output Structure

```
01_cleaning/
├── intermediate/   # SoupX outputs
├── output/         # Final cleaned H5AD per sample
└── qc/             # QC plots (PNG)

02_integration/
├── 03_output/      # Integrated H5AD
└── 04_output/      # Cell-typed H5AD, marker gene CSVs
```
