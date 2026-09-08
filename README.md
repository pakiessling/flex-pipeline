# flex-pipeline

A modular Snakemake pipeline for processing 10X Genomics FLEX single-cell RNA-seq data on SLURM HPC clusters.

## Prerequisites

| Software | Notes |
|---|---|
| Snakemake ≥ 8 | With `snakemake-executor-plugin-slurm` installed |
| Conda / Miniconda | For environment management |
| SLURM | For HPC job submission |

Install Snakemake and snakemake-executor-plugin-slurm in your PATH:
```bash
conda install -n base snakemake snakemake-executor-plugin-slurm
```
Make a copy of this pipeline:
```bash
git clone https://github.com/pakiessling/flex-pipeline
```
## Quick Start

### 1. Fill in your samples

Edit `samples.csv` (one row per sample). Two input modes are supported and can be mixed in the same run:

**CellRanger matrices** (runs SoupX + scDblFinder in step 1):
```csv
sample_id,raw_matrix_path,filtered_matrix_path
sample_A,/path/to/sampleA/raw_feature_bc_matrix.h5,/path/to/sampleA/filtered_feature_bc_matrix.h5
sample_B,/path/to/sampleB/raw_feature_bc_matrix.h5,/path/to/sampleB/filtered_feature_bc_matrix.h5
```

**Pre-processed h5ad files** (bypasses step 1, feeds directly into QC):
```csv
sample_id,h5ad_path
sample_A,/path/to/sampleA.h5ad
sample_B,/path/to/sampleB.h5ad
```

**Mixed** (both modes in one run):
```csv
sample_id,raw_matrix_path,filtered_matrix_path,h5ad_path
sample_A,/raw/A.h5,/filt/A.h5,
sample_B,,,/path/to/B.h5ad
```

When using pre-processed h5ad input, `scDblFinder.score`/`scDblFinder.class` annotations are optional — downstream steps handle their absence gracefully.

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
| 1 | `01_soupx_doublets.R` | Ambient RNA correction (SoupX) + doublet scoring (scDblFinder) | `soupx` |
| 2 | `02_qc.py` | QC annotations and technical sample eligibility | `qc` |
| 2 (optional) | `02_qc.py --diagnostics` | Per-sample clustering, markers, UMAP | `qc_diagnostics` |
| 3 | `03_integration.py` | Harmony batch correction, UMAP, PaCMAP | `integration` |
| 4 | `04_singleR.R` | Label transfer from reference dataset (SingleR) | `label_transfer` |
| 5 | `05_cytetype.py` | LLM-based cluster annotation (CyteType) | `cytetype` |
| 6 | `06_markers.py` | Marker genes via illico (Wilcoxon rank-sum) | `markers` |
| 7 | `07_report.py` | HTML summary report | `report` |

### Data flow

```
samples.csv
  ├─ CellRanger rows → [01_soupx_doublets]  results/intermediate/{sample}_cleaned.h5ad
  │                     (obs: scDblFinder.score, scDblFinder.class)
  └─ h5ad rows ──────→ copied directly  ──→ results/intermediate/{sample}_cleaned.h5ad
  → [02_qc]              results/per_sample/{sample}_clean.h5ad
                          (obs: cell_quality; counts retained)
                        + results/per_sample/{sample}_status.json
  → [sample_manifest]   results/integration/sample_inclusion.csv
  → [03_integration]     results/integration/integrated.h5ad
  → [04_singler]         results/annotation/integrated_labeled.h5ad   (optional)
  → [06_markers]         results/annotation/integrated_markers.h5ad
                        + results/markers/marker_genes_leiden_3.csv
  → [05_cytetype]        results/annotation/integrated_cytetype.h5ad  (optional)
                        + results/annotation/cytetype_annotation.json
  → [07_report]          results/report/report.html
```

### Doublet detection

Doublet detection is performed by [scDblFinder](https://bioconductor.org/packages/scDblFinder/) in step 1, on the SoupX-corrected counts before any QC filtering. Results are stored in two `.obs` columns that persist through the entire pipeline:

| Column | Type | Description |
|---|---|---|
| `scDblFinder.score` | float | Doublet score (higher = more likely doublet) |
| `scDblFinder.class` | string | `"singlet"` or `"doublet"` |

Doublets are **never hard-filtered** — they are carried as annotations so downstream analyses can choose how to handle them.

### Ambient RNA correction

Ambient RNA is corrected with [SoupX](github.com/constantAmateur/SoupX) 

This is only possible when unfiltered and filtered output of Cellranger is available.

The original counts are preserved in the object. 

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
  mad_threshold: 5                # MAD multiplier for per-sample QC outlier detection
  expected_doublet_rate: 0.05     # Passed as dbr to scDblFinder
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
| `r_soupx_doublets.yml` | Step 1 | R, SoupX, scDblFinder, anndataR |
| `sc.yml` | Steps 2, 3, 5, 6, 7 | scanpy, harmonypy, illico, cytetype, pacmap |
| `r_singler.yml` | Step 4 | R, SingleR, anndataR |

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
Set `steps.cytetype: true`. The markers step must also be enabled (`steps.markers: true`).
(Also see [CyteType docs](https://github.com/NygenAnalytics/CyteType)).

## Sample eligibility and cell retention

QC retains low-quality cells and scDblFinder doublets for filtering by the user
after the pipeline. `params.min_genes` now controls quality **annotation only**;
it no longer removes cells. MAD outliers are also annotations.

The only automatic cell removal in Python QC is a zero total count in the
selected expression matrix (after SoupX when available). Each removed original
barcode is recorded in `results/per_sample/{sample}_status.json`. Negative or
nonfinite expression values fail visibly as invalid inputs instead of being
silently removed. Preprocessed H5AD inputs should provide counts in
`after_soupx`, `b4_soupx`, `counts`, or otherwise `X`, in that priority order.

Samples are excluded from integration only when they cannot support the current
dimensional reduction: fewer than three nonzero-count cells, fewer than three
variable genes, or fewer than three variable genes after library-size
normalization and log transformation. These are technical feasibility checks,
not a biological sample-quality score. A biologically poor but technically
processable sample stays included.

Every assessed sample produces a QC H5AD and status JSON, including excluded
samples. The H5AD retains the full gene set and count layers for all remaining
cells. Integration reads an explicit manifest of the samples in the current
workflow, never all H5ADs found in a directory.

- `results/integration/sample_inclusion.csv` lists every requested sample,
  inclusion/exclusion reasons, and cell counts.
- The same table is embedded in `integrated.h5ad.uns["sample_inclusion"]` and
  shown in the final report.
- One eligible sample runs without Harmony; no eligible samples produces a clear
  integration error, leaving the manifest and QC results available.
- Known LOESS failures retry wider smoothing spans, then use Seurat HVGs on
  log-normalized data. The selected method and failures are recorded in
  `uns["hvg_selection"]`. Unknown errors still fail visibly.

Per-sample clustering, marker tables, and UMAP are now optional diagnostics,
separate from the QC files required by integration. Enable them with
`steps.qc_diagnostics: true` (default: false). They write
`results/qc/{sample}/{sample}_diagnostics.h5ad` and the existing plot/marker
paths. Ineligible samples receive a skipped diagnostic status. A diagnostic
failure cannot change sample eligibility; unexpected errors still fail that
diagnostic job. With Snakemake's keep-going setting, independent integration
jobs can continue.

After updating an existing run, rerun QC and integration to create the new
status files and manifest; existing SoupX outputs can be reused. The standard
`./launch.sh` command schedules missing status outputs automatically. No
environment-installation behavior is changed by this update.

Regression tests (in an environment with the Python pipeline dependencies and
pytest):

```bash
python -m pytest tests -q
```
