import argparse
import scanpy as sc
import celltypist
from celltypist import models
import random
import numpy as np
import os
from scipy import sparse

os.environ["PYTHONHASHSEED"] = "0"


np.random.seed(0)  # Numpy random
random.seed(0)  # Python random


def main(input_file, output_file):
    adata = sc.read_h5ad(input_file)
    # output_dir = os.path.dirname(output_file)

    # function for celltypist annotation
    def run_celltypist(adata, model_name, key_suffix):
        adata_celltypist = adata.copy()
        adata_celltypist.X = adata.layers["after_soupx"]
        sc.pp.normalize_per_cell(adata_celltypist, counts_per_cell_after=10**4)
        sc.pp.log1p(adata_celltypist)
        adata_celltypist.X = adata_celltypist.X.toarray()

        models.download_models(force_update=True)
        #  model = models.Model.load(model=model_name)

        predictions = celltypist.annotate(
            adata_celltypist, model=model_name, majority_voting=True
        )

        adata = predictions.to_adata()
        adata.obs[f"predicted_labels_{key_suffix}"] = adata.obs["predicted_labels"]
        adata.obs[f"majority_voting_{key_suffix}"] = adata.obs["majority_voting"]
        adata.obs[f"conf_score_{key_suffix}"] = adata.obs["conf_score"]

        sc.pl.umap(
            adata,
            color=[f"predicted_labels_{key_suffix}"],
            save=f"_predicted_labels_{key_suffix}.png",
        )
        sc.pl.umap(
            adata,
            color=[f"majority_voting_{key_suffix}"],
            save=f"_majority_voting_{key_suffix}.png",
        )

        return adata

    # Step 1: Immune cells model
    adata = run_celltypist(adata, "Immune_All_Low.pkl", "immune")

    # Step 2: Heart cells model
    adata = run_celltypist(adata, "Healthy_Adult_Heart.pkl", "heart")

    # Cell cycle analysis
    s_genes = [
        "MCM5",
        "PCNA",
        "TYMS",
        "FEN1",
        "MCM2",
        "MCM4",
        "RRM1",
        "UNG",
        "GINS2",
        "MCM6",
        "CDCA7",
        "DTL",
        "PRIM1",
        "UHRF1",
        "MLF1IP",
        "HELLS",
        "RFC2",
        "RPA2",
        "NASP",
        "RAD51AP1",
        "GMNN",
        "WDR76",
        "SLBP",
        "CCNE2",
        "UBR7",
        "POLD3",
        "MSH2",
        "ATAD2",
        "RAD51",
        "RRM2",
        "CDC45",
        "CDC6",
        "EXO1",
        "TIPIN",
        "DSCC1",
        "BLM",
        "CASP8AP2",
        "USP1",
        "CLSPN",
        "POLA1",
        "CHAF1B",
        "BRIP1",
        "E2F8",
    ]
    g2m_genes = [
        "HMGB2",
        "CDK1",
        "NUSAP1",
        "UBE2C",
        "BIRC5",
        "TPX2",
        "TOP2A",
        "NDC80",
        "CKS2",
        "NUF2",
        "CKS1B",
        "MKI67",
        "TMPO",
        "CENPF",
        "TACC3",
        "FAM64A",
        "SMC4",
        "CCNB2",
        "CKAP2L",
        "CKAP2",
        "AURKB",
        "BUB1",
        "KIF11",
        "ANP32E",
        "TUBB4B",
        "GTSE1",
        "KIF20B",
        "HJURP",
        "CDCA3",
        "HN1",
        "CDC20",
        "TTK",
        "CDC25C",
        "KIF2C",
        "RANGAP1",
        "NCAPD2",
        "DLGAP5",
        "CDCA2",
        "CDCA8",
        "ECT2",
        "KIF23",
        "HMMR",
        "AURKA",
        "PSRC1",
        "ANLN",
        "LBR",
        "CKAP5",
        "CENPE",
        "CTCF",
        "NEK2",
        "G2E3",
        "GAS2L3",
        "CBX5",
        "CENPA",
    ]

    sc.tl.score_genes_cell_cycle(adata, s_genes, g2m_genes)

    sc.pl.violin(
        adata,
        keys=["S_score", "G2M_score"],
        stripplot=False,
        save="_cell_cycle_score.png",
    )
    sc.pl.umap(adata, color=["phase"], save="_phase.png")
    
    sparse_X = sparse.csr_matrix(adata.X)
    adata.X = sparse_X

    # Save the final annotated data
    adata.write(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run cell typing on integrated scRNA-seq data"
    )
    parser.add_argument("--input_file", type=str, help="Path to the input h5ad file")
    parser.add_argument(
        "--output_file", type=str, help="Path to save the output h5ad file"
    )
    args = parser.parse_args()

    main(args.input_file, args.output_file)
