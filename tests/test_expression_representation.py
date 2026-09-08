from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from illico import asymptotic_wilcoxon

from test_sample_eligibility import counts, load_module
from processing_utils import run_pca_neighbors


def test_scaled_pca_does_not_change_expression():
    a = counts(np.random.default_rng(8).poisson(2, (20, 30)))
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    before = a.X.copy()
    a.var["highly_variable"] = [True] * 20 + [False] * 10
    run_pca_neighbors(a, scale=True)
    np.testing.assert_array_equal(a.X.toarray(), before.toarray())
    expected = a[:, a.var["highly_variable"]].copy()
    sc.pp.scale(expected, zero_center=False)
    sc.pp.pca(expected, n_comps=19)
    np.testing.assert_allclose(abs(a.obsm["X_pca"]), abs(expected.obsm["X_pca"]), atol=1e-4)
    assert a.varm["PCs"].shape == (30, 19)


def test_markers_read_logcounts_even_when_x_is_scaled(tmp_path):
    markers = load_module("markers", "06_markers.py")
    a = counts([[2, 10], [3, 12], [4, 9], [20, 10], [25, 11], [30, 12]])
    a.obs["group"] = pd.Categorical(["A"] * 3 + ["B"] * 3)
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    a.layers["logcounts"] = a.X.copy()
    expected = asymptotic_wilcoxon(a, group_keys="group", reference=None, is_log1p=True)
    sc.pp.scale(a, zero_center=False)
    source = tmp_path / "input.h5ad"
    a.write_h5ad(source)
    args = SimpleNamespace(input_file=str(source), output_file=str(tmp_path / "output.h5ad"),
                           expression_layer="logcounts", group_key="group", is_log1p=True,
                           markers_dir=str(tmp_path), n_top=2, pval_cutoff=1.0)
    markers.main(args)
    actual = pd.read_csv(tmp_path / "all_genes_group.csv")
    np.testing.assert_allclose(actual["fold_change"], expected["fold_change"])
    output = ad.read_h5ad(args.output_file)
    np.testing.assert_array_equal(output.X.toarray(), a.X.toarray())
    assert output.uns["rank_genes_groups"]["params"]["layer"] == "logcounts"
