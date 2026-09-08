"""Optional PaCMAP and standalone plots, independent of the integration output."""
import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pacmap
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def main(args):
    if Path(args.input).resolve() == Path(args.output).resolve():
        raise ValueError("Diagnostics must not overwrite the core integrated object")
    adata = sc.read_h5ad(args.input)
    leiden_resolutions = [float(r) for r in args.leiden_resolutions.split()]
    plot_dir = str(Path(args.output).parent)
    Path(plot_dir).mkdir(parents=True, exist_ok=True)
    np.random.seed(0)
    logger.info("Computing PaCMAP …")
    if adata.n_obs >= 20:
        embedding = pacmap.PaCMAP(n_neighbors=min(10, adata.n_obs - 1))
        adata.obsm["X_pacmap"] = embedding.fit_transform(adata.obsm["X_pca_harmony"])
        adata.uns["pacmap_status"] = "completed"
    else:
        adata.uns["pacmap_status"] = "skipped: fewer than 20 cells"
        logger.info(adata.uns["pacmap_status"])

    logger.info("Saving UMAP plots …")
    sc.settings.figdir = plot_dir

    ri = np.random.permutation(adata.n_obs)
    sc.pl.umap(adata[ri, :], color="Sample", save="_by_sample.png", show=False)

    cluster_cols = [f"leiden_{str(r).replace('.', '_')}" for r in leiden_resolutions]
    sc.pl.umap(
        adata,
        color=cluster_cols,
        legend_loc="on data",
        save="_clusters.png",
        show=False,
    )

    qc_cols = [c for c in ["scDblFinder.class", "cell_quality"] if c in adata.obs.columns]
    if qc_cols:
        sc.pl.umap(adata[ri, :], color=qc_cols, save="_qc.png", show=False)

    adata.uns["integration_diagnostics"] = {
        "completed": True, "pacmap_status": adata.uns["pacmap_status"],
        "pacmap_version": pacmap.__version__,
    }
    adata.write_h5ad(args.output)
    plt.close("all")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--leiden_resolutions", default="1.5 3.0")
    main(parser.parse_args())
