"""Technical feasibility checks; these are not biological quality filters."""

import logging

import numpy as np
import scanpy as sc
from scipy import sparse

logger = logging.getLogger(__name__)


class InsufficientData(ValueError):
    """The available data cannot support dimensional reduction."""


def variable_mask(x):
    if x.shape[0] < 2:
        return np.zeros(x.shape[1], dtype=bool)
    x = x.astype(np.float64)
    mean = np.asarray(x.mean(axis=0)).ravel()
    second = np.asarray(x.multiply(x).mean(axis=0)).ravel() if sparse.issparse(x) else np.mean(x * x, axis=0)
    return (second - mean * mean) > np.maximum(second, 1) * 1e-12


def prepare_counts(adata, sample):
    """Preserve counts and remove only zero-total cells; invalid values are errors."""
    source = next((k for k in ("after_soupx", "b4_soupx", "counts") if k in adata.layers), None)
    if source:
        adata.X = adata.layers[source].copy()
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"[{sample}] Expected finite, nonnegative counts; check the input representation")
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    removed = adata.obs_names[totals == 0].tolist()
    adata = adata[totals > 0].copy()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    n_variable = int(variable_mask(adata.X).sum())
    # Counts differing only in library size become constant after normalization.
    normalized = adata.X.astype(np.float64).copy()
    if adata.n_obs:
        nonzero_totals = totals[totals > 0]
        factors = np.median(nonzero_totals) / nonzero_totals
        if sparse.issparse(normalized):
            normalized = normalized.multiply(factors[:, None]).tocsr()
            normalized.data = np.log1p(normalized.data)
        else:
            normalized = np.log1p(normalized * factors[:, None])
    n_variable_normalized = int(variable_mask(normalized).sum())
    reasons = []
    if adata.n_obs < 3:
        reasons.append("fewer than 3 nonzero-count cells")
    if n_variable < 3:
        reasons.append("fewer than 3 variable genes")
    if n_variable >= 3 and n_variable_normalized < 3:
        reasons.append("fewer than 3 variable genes after normalization")
    status = {
        "sample": sample,
        "integration_eligible": not reasons,
        "reason": "; ".join(reasons),
        "n_cells_input": adata.n_obs + len(removed),
        "n_cells_output": adata.n_obs,
        "n_zero_count_cells_removed": len(removed),
        "removed_cell_ids": removed,
        "n_variable_genes": n_variable,
        "n_variable_genes_normalized": n_variable_normalized,
    }
    return adata, status


def is_loess_error(error):
    message = str(error).lower()
    return any(s in message for s in (
        "reciprocal condition number", "singular", "extrapolation",
        "svddc", "span is too small", "near neighborhood", "zero-width neighborhood",
    ))


def select_hvgs(adata, n_top_genes=4000, batch_key=None):
    """Select on a copy, leaving counts intact, with logged LOESS fallbacks."""
    informative = variable_mask(adata.X)
    if adata.n_obs < 3 or informative.sum() < 3:
        raise InsufficientData("Need at least 3 cells and 3 variable genes")
    work = adata[:, informative].copy()
    n_top = min(n_top_genes, work.n_vars)
    if n_top < 3:
        raise ValueError("n_top_genes must be at least 3")
    failures = []
    flavor = "seurat_v3_paper" if batch_key else "seurat_v3"
    for span in (0.3, 0.5, 1.0):
        try:
            sc.pp.highly_variable_genes(work, flavor=flavor, span=span,
                                       n_top_genes=n_top, batch_key=batch_key)
            method = f"{flavor}; span={span}"
            break
        except ValueError as error:
            if not is_loess_error(error):
                raise
            failures.append(str(error))
            logger.warning("HVG %s span=%s failed: %s", flavor, span, error)
    else:
        sc.pp.normalize_total(work)
        sc.pp.log1p(work)
        sc.pp.highly_variable_genes(work, flavor="seurat", n_top_genes=n_top,
                                   batch_key=batch_key)
        method = "seurat (log-normalized counts)"
    mask = np.zeros(adata.n_vars, dtype=bool)
    mask[informative] = work.var["highly_variable"].to_numpy()
    if mask.sum() < 3:
        raise InsufficientData("HVG selection yielded fewer than 3 genes")
    adata.var["highly_variable"] = mask
    adata.uns["hvg_selection"] = {"method": method, "failures": failures}
    logger.info("HVG selection: %s (%d genes)", method, mask.sum())


def run_pca_neighbors(adata):
    n_comps = min(50, adata.n_obs - 1, int(adata.var["highly_variable"].sum()) - 1)
    if n_comps < 2:
        raise InsufficientData("Fewer than 2 feasible principal components")
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1), use_rep="X_pca")
