# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a **Snakemake-based bioinformatics pipeline** for processing 10X Genomics FLEX single-cell RNA-seq data. The pipeline is unified and fully modular — all steps run from CellRanger output to final report, controlled by a single config file.

## Running the Pipeline

**Start pipeline (SLURM HPC):**
```bash
./launch.sh
```

**Dry-run (preview without submitting):**
```bash
./launch.sh --dry-run
```

**Run locally (testing):**
```bash
./launch.sh --local
```

## Architecture

### Data Flow

```
samples.csv + CellRanger H5 matrices
  → [01_soupx.py]      Ambient RNA correction via R/SoupX (optional)
  → [02_qc.py]         QC metrics, doublet detection, clustering per sample
  → [03_integration.py] Merge samples, Harmony batch correction, UMAP/PaCMAP
  → [04_singleR.R]     Label transfer from reference h5ad via SingleR (optional)
  → [06_markers.py]    Wilcoxon rank-sum marker genes via illico
  → [05_cytetype.py]   LLM-based cluster annotation via CyteType (optional)
  → [07_report.py]     HTML summary report
```

### Key Design Decisions

- **Modular steps**: each step is toggled in `config/config.yaml` under `steps:`. Disabled steps either produce passthrough outputs or are excluded from the DAG.
- **Single config**: all analysis parameters, step toggles, and cluster settings live in `config/config.yaml` and `config/cluster.yaml`. No hardcoded values in scripts.
- **Cells marked not removed** during QC (`02_qc.py`). Hard filter for doublets happens in `03_integration.py`.
- **SoupX graceful fallback**: `01_soupx.py` falls back to the unmodified filtered matrix on failure.
- **Reproducibility**: all scripts set `np.random.seed(0)`, `random.seed(0)`, `PYTHONHASHSEED=0`, and write software versions + parameters to `adata.uns["pipeline_log"]`.
- **AnnData layers**: `b4_soupx` and `after_soupx`. Integration checks which layer to use.
- **illico replaces R/Presto**: marker genes are computed via `illico.asymptotic_wilcoxon` and stored in scanpy's `rank_genes_groups` format so CyteType can consume them.
- **Cell cycle genes** are loaded from `config/cell_cycle_genes.json` (not hardcoded in scripts).

### Directory Structure

```
flex-pipeline/
├── config/
│   ├── config.yaml          # Master config: steps, params, output dirs
│   ├── cluster.yaml         # SLURM account, partition, conda path
│   └── cell_cycle_genes.json
├── workflow/
│   ├── Snakefile            # Unified pipeline DAG
│   ├── scripts/
│   │   ├── 01_soupx.py
│   │   ├── 02_qc.py
│   │   ├── 03_integration.py
│   │   ├── 04_singleR.R
│   │   ├── 05_cytetype.py
│   │   ├── 06_markers.py
│   │   └── 07_report.py
│   └── environments/
│       ├── py_r.yml         # Step 1: Python + R + SoupX
│       ├── sc.yml           # Steps 2,3,5,6,7: scanpy + illico + cytetype
│       └── r_singler.yml    # Step 4: R + SingleR
├── profiles/
│   └── config.yaml          # snakemake-executor-plugin-slurm settings
├── samples.csv              # Fill in paths to CellRanger outputs
├── launch.sh                # Entry point
├── jobscript.sh             # SLURM job template
└── results/                 # Created at runtime
    ├── intermediate/        # Per-sample SoupX outputs
    ├── per_sample/          # Per-sample QC h5ad files
    ├── qc/                  # QC plots
    ├── integration/         # Integrated h5ad
    ├── annotation/          # SingleR + CyteType outputs
    ├── markers/             # Marker gene CSVs
    └── report/              # HTML report
```

### Conda Environments

| File | Used by | Key packages |
|------|---------|--------------|
| `workflow/environments/py_r.yml` | `run_soupx` / `skip_soupx` | Python 3.11, R 4.3, SoupX, rpy2, anndata2ri |
| `workflow/environments/sc.yml` | `run_qc`, `run_integration`, `run_markers`, `run_cytetype`, `run_report` | scanpy, harmonypy, illico, cytetype, pacmap |
| `workflow/environments/r_singler.yml` | `run_label_transfer` | R 4.3, SingleR, anndataR |

### SLURM Configuration

Uses `snakemake-executor-plugin-slurm` (configured in `profiles/config.yaml`).
Default resources: 50 GB RAM, 1440 min walltime, 4 CPUs.
Integration and markers rules request 180 GB / 12 CPUs.
Account and partition are read from `config/cluster.yaml` by `launch.sh`.

### Input Format

`samples.csv` columns:
- `sample_id` — unique identifier
- `raw_matrix_path` — path to CellRanger raw matrix H5 file
- `filtered_matrix_path` — path to CellRanger filtered matrix H5 file
