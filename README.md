# flex-pipeline

A modular Snakemake pipeline for processing 10X Genomics FLEX single-cell RNA-seq data on SLURM HPC clusters.

## Prerequisites

| Software | Notes |
|---|---|
| Snakemake ≥ 8 | With `snakemake-executor-plugin-slurm` installed |
| Conda / Miniconda | For environment management |
| SLURM | For HPC job submission |

## Quick Start

### 1. Fill in your samples

Edit `samples.csv` (one row per sample):

```csv
sample_id,raw_matrix_path,filtered_matrix_path
sample_A,/path/to/sampleA/raw_feature_bc_matrix.h5,/path/to/sampleA/filtered_feature_bc_matrix.h5
sample_B,/path/to/sampleB/raw_feature_bc_matrix.h5,/path/to/sampleB/filtered_feature_bc_matrix.h5
```

### 2. Configure the pipeline

Edit **`config/config.yaml`** to set which steps to run and tweak analysis parameters.

Edit **`config/cluster.yaml`** to set your SLURM account, partition, and conda path.

### 3. Run

```bash
# Preview jobs without submitting (dry-run)
./launch.sh --dry-run

# Submit to SLURM
./launch.sh

# Run locally (no SLURM, useful for testing)
./launch.sh --local
```

Logs are written to `logs/pipeline_YYYYMMDD_HHMMSS.log` and to `logs/<step>/` per rule.

---

## Pipeline Steps

All steps can be toggled on/off in `config/config.yaml` under `steps:`.

| Step | Script | Description | Toggle |
|---|---|---|---|
| 1 | `01_soupx.py` | Ambient RNA removal (SoupX) | `soupx` |
| 2 | `02_qc.py` | Per-sample QC, doublet detection, clustering | `qc` |
| 3 | `03_integration.py` | Harmony batch correction, UMAP, PaCMAP | `integration` |
| 4 | `04_singleR.R` | Label transfer from reference dataset (SingleR) | `label_transfer` |
| 5 | `05_cytetype.py` | LLM-based cluster annotation (CyteType) | `cytetype` |
| 6 | `06_markers.py` | Marker genes via illico (Wilcoxon rank-sum) | `markers` |
| 7 | `07_report.py` | HTML summary report | `report` |

### Data flow

```
samples.csv + CellRanger H5 matrices
  → [01_soupx]   results/intermediate/{sample}_cleaned.h5ad
  → [02_qc]      results/per_sample/{sample}_clean.h5ad
  → [03_intgr]   results/integration/integrated.h5ad
  → [04_singler] results/annotation/integrated_labeled.h5ad   (optional)
  → [06_markers] results/annotation/integrated_markers.h5ad
               + results/markers/marker_genes_leiden_3.csv
  → [05_cytetype] results/annotation/integrated_cytetype.h5ad (optional)
                + results/annotation/cytetype_annotation.json
  → [07_report]  results/report/report.html
```

---

## Configuration Reference

### `config/config.yaml`

```yaml
samples: "samples.csv"

steps:
  soupx: true           # Skip if you don't have raw matrices
  label_transfer: false # Requires reference_h5ad to be set
  cytetype: false       # Requires CyteType API access

params:
  n_top_genes: 4000
  leiden_resolutions: [1.5, 3.0]
  reference_h5ad: "/path/to/reference.h5ad"  # For label_transfer step
  cytetype_study_context: "Human heart tissue, 10X FLEX ..."
```

### `config/cluster.yaml`

```yaml
slurm_account:    "your_account"
slurm_partition:  "your_partition"
conda_root:       ""   # Leave blank to auto-detect $HOME/miniconda3
snakemake_env:    ""   # Leave blank if snakemake is already in PATH
```

---

## Environments

Three conda environments are used (defined in `workflow/environments/`):

| File | Used by | Key packages |
|---|---|---|
| `py_r.yml` | Step 1 (SoupX) | Python + R + rpy2 + SoupX |
| `sc.yml` | Steps 2, 3, 5, 6, 7 | scanpy, harmonypy, illico, cytetype |
| `r_singler.yml` | Step 4 (SingleR) | R + SingleR + anndataR |

Build environments with:
```bash
conda env create -f workflow/environments/sc.yml
```

---

## Optional steps — setup notes

### Label transfer (SingleR)
Set `steps.label_transfer: true` and provide a reference h5ad:
```yaml
params:
  reference_h5ad: "/path/to/reference.h5ad"
  reference_label_column: "cell_type"
```

### CyteType annotation
Set `steps.cytetype: true`. CyteType requires an API key (see [CyteType docs](https://github.com/NygenAnalytics/CyteType)).
The markers step must also be enabled (`steps.markers: true`).

---

## Notes on the old two-phase structure

The previous `01_cleaning/` and `02_integration/` directories are retained for reference but superseded by this unified pipeline. They can be safely deleted once you have validated the new pipeline on your data.
