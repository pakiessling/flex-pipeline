import argparse
import logging
import os

import anndata2ri
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
import scanpy as sc

# Set up logging
rcb.logger.setLevel(logging.INFO)
ro.pandas2ri.activate()
anndata2ri.activate()

# Argument parser
parser = argparse.ArgumentParser(description="Run SoupX for ambient RNA removal.")
parser.add_argument("--sample", required=True, help="Sample identifier")
parser.add_argument("--output", required=True, help="Output path for cleaned data")
parser.add_argument(
    "--raw_matrix", required=True, help="Path to the raw feature matrix"
)
parser.add_argument(
    "--filtered_matrix", required=True, help="Path to the filtered feature matrix"
)
args = parser.parse_args()

# Load the filtered matrix
adata = sc.read_10x_h5(args.filtered_matrix)
adata.var_names_make_unique()
adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp)
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added="soupx_groups")

# Preprocess variables for SoupX
soupx_groups = adata_pp.obs["soupx_groups"]
del adata_pp

# Load the raw matrix and filter only common genes
adata_raw = sc.read_10x_h5(args.raw_matrix)
adata_raw.var_names_make_unique()
common_genes = adata_raw.var_names.intersection(adata.var_names)
common_genes = common_genes.sort_values()
adata = adata[:, common_genes].copy()
adata_raw = adata_raw[:, common_genes].copy()
data_tod = adata_raw.X.T
del adata_raw

cells = adata.obs_names
genes = adata.var_names
data = adata.X.T

# Load Python variables into R (for SoupX processing)
ro.globalenv["data"] = data
ro.globalenv["data_tod"] = data_tod
ro.globalenv["genes"] = genes
ro.globalenv["cells"] = cells
ro.globalenv["soupx_groups"] = soupx_groups

# Initialize a flag to check if SoupX succeeded
soupx_success = True

try:
    # R code execution with rpy2
    ro.r("""
       library(SoupX)
       library(Matrix)

       # specify row and column names of data
       rownames(data) <- genes
       colnames(data) <- cells

       # ensure correct sparse format for table of counts and table of droplets
       data <- as(data, "sparseMatrix")
       data_tod <- as(data_tod, "sparseMatrix")

       # Generate SoupChannel Object for SoupX
       sc <- SoupChannel(data_tod, data, calcSoupProfile = FALSE)

       # Add extra meta data to the SoupChannel object
       soupProf <- data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
       sc <- setSoupProfile(sc, soupProf)

       # Set cluster information in SoupChannel
       sc <- setClusters(sc, soupx_groups)

       # Estimate contamination fraction
       sc <- autoEstCont(sc, doPlot=FALSE)

       # Infer corrected table of counts and round to integer
       out <- adjustCounts(sc, roundToInt = TRUE)
    """)

    # Retrieve the adjusted counts from R
    out = ro.globalenv["out"]

    # Store SoupX corrected counts in AnnData
    adata.layers["b4_soupx"] = adata.X.copy()
    adata.layers["after_soupx"] = out.T
    adata.X = adata.layers["after_soupx"].copy()

except Exception as e:
    adata.layers["b4_soupx"] = adata.X.copy()
    soupx_success = False
    print(f"SoupX processing failed: {e}")
    logging.error(f"SoupX processing failed: {e}")

# Print total number of genes before and after filtering
print(f"Total number of genes: {adata.n_vars}")

os.makedirs(os.path.dirname(args.output), exist_ok=True)
# Write the cleaned AnnData object to the specified output path
adata.write(args.output)
print(
    "SoupX processing completed successfully"
    if soupx_success
    else "SoupX was skipped, data saved without it."
)

# Add a summary of the changes
if soupx_success:
    before_sum = adata.layers["b4_soupx"].sum()
    after_sum = adata.layers["after_soupx"].sum()
    percent_change = (after_sum - before_sum) / before_sum * 100
    print(f"Total count before SoupX: {before_sum}")
    print(f"Total count after SoupX: {after_sum}")
