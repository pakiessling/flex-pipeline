import argparse
import os
import warnings

import doubletdetection
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation as mad

os.environ["PYTHONHASHSEED"] = "0"
import random

np.random.seed(0)  # Numpy random
random.seed(0)  # Python random

# mute noise during marker gene calculation
warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

# Argument parser
parser = argparse.ArgumentParser(
    description="Quality control for single-cell RNA-seq data."
)
parser.add_argument("--input", required=True, help="Path to the input .h5ad file")
parser.add_argument(
    "--sample", required=True, help="Sample identifier for naming outputs"
)
parser.add_argument("--output", required=True, help="Output path for cleaned data")
parser.add_argument(
    "--qc_folder", required=True, help="Output path for marker genes and UMAPs"
)
args = parser.parse_args()

# Load anndata and make a new obs column with sample metadata
adata = sc.read_h5ad(args.input)
adata.obs["Sample"] = args.sample
adata.obs.index = adata.obs.index + "_" + adata.obs.Sample

print("adata loaded")


def qc(adata):
    adata.var["MT"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["MT"], inplace=True, percent_top=[20], log1p=True
    )
    return adata


adata = qc(adata)


# Function to remove outlier cells based on median absolute deviations
def mad_outlier(adata, metric, nmads, upper_only=False):
    M = adata.obs[metric]
    if not upper_only:
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))
    return M > np.median(M) + nmads * mad(M)


# Removes outlier cells from an individual anndata object
def pp(adata):
    sc.pp.filter_cells(adata, min_genes=2)
    bool_vector = (
        mad_outlier(adata, "log1p_total_counts", 5)
        + mad_outlier(adata, "log1p_n_genes_by_counts", 5)
        + mad_outlier(adata, "pct_counts_in_top_20_genes", 5)
        # + mad_outlier(adata, "pct_counts_MT", 3, upper_only=True) # MT RNA is not removed for this dataset
    )

    # Add a new column in adata.obs to label low-quality cells
    adata.obs["cell_quality"] = "high-quality"  # Default to "high-quality"
    adata.obs.loc[bool_vector, "cell_quality"] = "low-quality"  # Mark low-quality cells

    # Optionally filter out low-quality cells
    # adata = adata[~bool_vector]

    return adata


adata = pp(adata)

print(
    adata.obs["cell_quality"].value_counts()
)  # to check the amount of low-quality cells

print("QC done")

# Doublet detection
clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)


def dd(adata):
    sc.pp.scrublet(adata, expected_doublet_rate=0.2)
    doublets = clf.fit(adata.X).predict(p_thresh=1e-10, voter_thresh=0.5)
    doublet_score = clf.doublet_score()
    adata.obs["clf_doublet"] = doublets
    adata.obs["clf_score"] = doublet_score
    return adata


adata = dd(adata)

print("Doublet detection done")

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=1.5, key_added="leiden_1_5")
sc.tl.leiden(adata, resolution=3, key_added="leiden_3")

sc.tl.rank_genes_groups(
    adata,
    "leiden_1_5",
    method="wilcoxon",
    tie_correct=True,
    key_added="rank_genes_groups_1_5",
)
sc.tl.rank_genes_groups(
    adata,
    "leiden_3",
    method="wilcoxon",
    tie_correct=True,
    key_added="rank_genes_groups_3",
)

# need to do this after rank_genes_groups as it can introduce negative values
sc.pp.scale(adata, zero_center=False)
sc.pp.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata)

sc.tl.umap(adata)

print("Creating UMAPs ...")
os.makedirs(f"{args.qc_folder}/{args.sample}", exist_ok=True)


sc.pl.umap(
    adata,
    color=[
        "predicted_doublet",
        "clf_doublet",
        "leiden_1_5",
        "leiden_3",
        "log1p_total_counts",
        "log1p_n_genes_by_counts",
        "pct_counts_MT",
    ],
    legend_loc="on data",
    return_fig=True,
    vmax="p99",
    sort_order=False,
)
plt.savefig(
    f"{args.qc_folder}/{args.sample}/{args.sample}_umap.png",
    bbox_inches="tight",
    dpi=200,
)
plt.close()


top_markers = pd.DataFrame(adata.uns["rank_genes_groups_1_5"]["names"]).head(100)
top_markers.to_csv(f"{args.qc_folder}/{args.sample}/{args.sample}_marker_1_5.csv")

top_markers_2 = pd.DataFrame(adata.uns["rank_genes_groups_3"]["names"]).head(100)
top_markers_2.to_csv(f"{args.qc_folder}/{args.sample}/{args.sample}_marker_3.csv")

# Write final result
os.makedirs(os.path.dirname(args.output), exist_ok=True)
adata.write_h5ad(args.output)

print("All done")
