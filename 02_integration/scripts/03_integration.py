import argparse
import os

import numpy as np
import pacmap
import scanpy as sc

os.environ["PYTHONHASHSEED"] = "0"
import random

np.random.seed(0)  # Numpy random
random.seed(0)  # Python random


def main(input_dir, output_file):
    # Step 1: Load and concat all files
    file_paths = [os.path.join(input_dir, x) for x in os.listdir(input_dir) if x.endswith(".h5ad")]

    def load(path):
        return sc.read_h5ad(path)

    adatas = [load(x) for x in file_paths]

    for adata in adatas:
        sample_name = adata.obs["Sample"].unique()
        # Check if SoupX was run by checking for the "after_soupx" layer
        if "after_soupx" in adata.layers:
            adata.X = adata.layers["after_soupx"].copy()
            adata.obs["SoupX_run"] = True
            print(f"Used 'after_soupx' for sample: {sample_name}")
        elif "b4_soupx" in adata.layers:
            adata.X = adata.layers["b4_soupx"].copy()
            adata.obs["SoupX_run"] = False
            print(f"Used 'b4_soupx' for sample: {sample_name}")
        else:
            adata.obs["SoupX_run"] = False
            print(f"Neither 'after_soupx' nor 'b4_soupx' layer found for sample: {sample_name}")

    adata = sc.concat(adatas, join="outer", merge="same")

    # remove trash after annotation
    # clusters_to_remove = ["Trash", "TRash", "trash"]
    # adata = adata[~adata.obs["cell_type"].isin(clusters_to_remove)].copy()
    adata = adata[adata.obs.predicted_doublet != True].copy()
    # Remove doublets detected by scrublet
    # adata = adata[~adata.obs.predicted_doublet].copy()

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3_paper",
        batch_key="Sample",
    )
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, zero_center=False)
    sc.pp.pca(adata)
    sc.external.pp.harmony_integrate(adata, max_iter_harmony=100, key=["Sample"])
    sc.pp.neighbors(adata, use_rep="X_pca_harmony")

    sc.tl.leiden(adata, resolution=1.5, key_added="leiden_1_5")
    sc.tl.leiden(adata, resolution=3, key_added="leiden_3")

    sc.tl.umap(adata)

    embedding = pacmap.PaCMAP()
    adata.obsm["X_pacmap"] = embedding.fit_transform(adata.obsm["X_pca_harmony"])

    # Generate plots
    ri = np.random.permutation(list(range(adata.shape[0])))
    sc.pl.umap(adata[ri, :], color="Sample", save="_intgr_all_heart_samples.png")
    sc.pl.umap(
        adata,
        color=["leiden_1_5", "leiden_3"],
        legend_loc="on data",
        save="_clusters.png",
    )

    # Save the integrated data
    adata.write(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Integrate scRNA-seq data")
    parser.add_argument("--input_dir", type=str, help="Directory containing input h5ad files")
    parser.add_argument("--output_file", type=str, help="Path to save the integrated h5ad file")
    args = parser.parse_args()

    main(args.input_dir, args.output_file)
