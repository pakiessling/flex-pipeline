import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
from processing_utils import prepare_counts, select_hvgs
from sample_manifest import build_manifest


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


qc = load_module("qc", "02_qc.py")
integration = load_module("integration", "03_integration.py")


def counts(x):
    a = ad.AnnData(sparse.csr_matrix(x, dtype=np.float64))
    a.obs_names = [f"cell{i}" for i in range(a.n_obs)]
    a.var_names = [f"gene{i}" for i in range(a.n_vars)]
    return a


def run_qc(tmp_path, sample, a):
    source = tmp_path / f"{sample}.h5ad"
    output = tmp_path / f"{sample}_clean.h5ad"
    status = tmp_path / f"{sample}.json"
    a.write_h5ad(source)
    qc.main(SimpleNamespace(
        input=str(source), output=str(output), sample=sample, diagnostics=False,
        status_output=str(status), qc_folder=str(tmp_path / "plots"),
        min_genes=2, mad_threshold=5, leiden_resolutions="1.0",
    ))
    return output, status


def test_preserves_low_quality_doublets_and_counts(tmp_path):
    a = counts([[0, 0, 0, 0], [1, 0, 0, 0], [1, 2, 0, 3], [2, 1, 4, 0]])
    a.obs["scDblFinder.class"] = ["singlet", "doublet", "singlet", "singlet"]
    output, status = run_qc(tmp_path, "A", a)
    result = ad.read_h5ad(output)
    audit = json.loads(status.read_text())
    assert audit["integration_eligible"]
    assert audit["removed_cell_ids"] == ["cell0"]
    assert result.n_obs == 3
    assert result.obs.loc["cell1_A", "cell_quality"] == "low-quality"
    assert result.obs.loc["cell1_A", "scDblFinder.class"] == "doublet"
    np.testing.assert_array_equal(result.layers["counts"].toarray(), a.X.toarray()[1:])
    np.testing.assert_array_equal(result.X.toarray(), a.X.toarray()[1:])


@pytest.mark.parametrize("x", [np.zeros((4, 5)), np.ones((4, 5)), np.ones((2, 5)), np.zeros((3, 0))])
def test_unusable_samples_write_valid_outputs(tmp_path, x):
    output, status = run_qc(tmp_path, "bad", counts(x))
    assert not json.loads(status.read_text())["integration_eligible"]
    assert not ad.read_h5ad(output).uns["sample_status"]["integration_eligible"]


@pytest.mark.parametrize("value", [np.nan, np.inf, -1])
def test_invalid_input_fails_visibly(value):
    a = counts([[value, 1], [2, 1]])
    with pytest.raises(ValueError, match="finite, nonnegative"):
        prepare_counts(a, "invalid")


def test_hvg_fallback_keeps_original_counts(monkeypatch):
    a = counts(np.random.default_rng(4).poisson(2, (30, 80)))
    original = a.X.copy()
    real = qc.sc.pp.highly_variable_genes
    attempts = []

    def fail_loess(work, **kwargs):
        attempts.append(kwargs.copy())
        if kwargs["flavor"].startswith("seurat_v3"):
            raise ValueError(b"reciprocal condition number 1.7043e-14")
        assert "log1p" in work.uns
        return real(work, **kwargs)

    monkeypatch.setattr(qc.sc.pp, "highly_variable_genes", fail_loess)
    a.obs["Sample"] = ["A"] * 15 + ["B"] * 15
    select_hvgs(a, n_top_genes=20, batch_key="Sample")
    assert [kw.get("span") for kw in attempts[:3]] == [0.3, 0.5, 1.0]
    assert all(kw["batch_key"] == "Sample" for kw in attempts)
    assert a.uns["hvg_selection"]["method"].startswith("seurat (")
    np.testing.assert_array_equal(a.X.toarray(), original.toarray())


def test_unknown_hvg_error_is_not_swallowed(monkeypatch):
    a = counts(np.random.default_rng(2).poisson(2, (10, 20)))

    def fail(*args, **kwargs):
        raise ValueError("unexpected bug")

    monkeypatch.setattr(qc.sc.pp, "highly_variable_genes", fail)
    with pytest.raises(ValueError, match="unexpected bug"):
        select_hvgs(a)


def test_manifest_excludes_bad_and_ignores_stale_files(tmp_path):
    good = run_qc(tmp_path, "good", counts(np.random.default_rng(2).poisson(2, (10, 30))))
    bad = run_qc(tmp_path, "bad", counts(np.zeros((4, 30))))
    counts(np.ones((5, 30))).write_h5ad(tmp_path / "stale.h5ad")
    manifest = tmp_path / "sample_inclusion.csv"
    build_manifest([good[0], bad[0]], [good[1], bad[1]], manifest)
    rows = pd.read_csv(manifest, keep_default_na=False)
    loaded = integration.load_samples(rows)
    assert len(loaded) == 1
    assert loaded[0].obs["Sample"].unique().tolist() == ["good"]
    assert rows["status"].tolist() == ["included", "excluded"]
    assert rows["n_cells_integrated"].tolist() == [10, 0]


def test_no_eligible_samples_keeps_manifest(tmp_path):
    bad = run_qc(tmp_path, "bad", counts(np.zeros((3, 5))))
    manifest = tmp_path / "sample_inclusion.csv"
    build_manifest([bad[0]], [bad[1]], manifest)
    with pytest.raises(ValueError, match="No samples are eligible"):
        integration.load_samples(pd.read_csv(manifest, keep_default_na=False))
    assert manifest.exists()


@pytest.mark.parametrize("n_samples", [1, 2])
def test_small_integration_round_trip(tmp_path, n_samples):
    rng = np.random.default_rng(12)
    outputs = [run_qc(tmp_path, f"S{i}", counts(rng.poisson(2, (3, 30))))
               for i in range(n_samples)]
    bad = run_qc(tmp_path, "excluded", counts(np.zeros((3, 30))))
    outputs.append(bad)
    manifest = tmp_path / "sample_inclusion.csv"
    build_manifest([v[0] for v in outputs], [v[1] for v in outputs], manifest)
    output = tmp_path / "integrated.h5ad"
    integration.main(SimpleNamespace(
        manifest=str(manifest), output_file=str(output), n_top_genes=20,
        leiden_resolutions="1.0", max_iter_harmony=5, cell_cycle_genes="",
    ))
    result = ad.read_h5ad(output)
    assert result.n_obs == 3 * n_samples
    assert np.isfinite(result.obsm["X_umap"]).all()
    assert result.uns["pipeline_log"]["integration"]["harmony_run"] == (n_samples > 1)
    assert result.uns["sample_inclusion"]["status"].tolist() == ["included"] * n_samples + ["excluded"]
    assert result.layers["counts"].shape == result.shape


def test_skipped_diagnostics_do_not_block_qc(tmp_path, monkeypatch):
    output, _ = run_qc(tmp_path, "A", counts(np.random.default_rng(2).poisson(2, (8, 30))))

    def fail(*args):
        raise np.linalg.LinAlgError("test numerical failure")

    monkeypatch.setattr(qc, "cluster_and_embed", fail)
    diagnostic = tmp_path / "diagnostic.h5ad"
    qc.main(SimpleNamespace(
        input=str(output), output=str(diagnostic), sample="A", diagnostics=True,
        qc_folder=str(tmp_path), leiden_resolutions="1.0",
    ))
    assert ad.read_h5ad(output).uns["sample_status"]["integration_eligible"]
    assert ad.read_h5ad(diagnostic).uns["diagnostics"]["status"] == "skipped"


def test_library_size_only_variation_is_excluded():
    a = counts([[1, 2, 3], [2, 4, 6], [3, 6, 9]])
    result, status = prepare_counts(a, "proportional")
    assert result.n_obs == 3
    assert not status["integration_eligible"]
    assert "after normalization" in status["reason"]


def test_dense_count_layers_are_preserved():
    a = counts([[1, 0, 2], [0, 3, 1], [3, 2, 0]])
    a.X = a.X.toarray()
    a.layers["b4_soupx"] = a.X.copy()
    result, status = prepare_counts(a, "dense")
    assert status["integration_eligible"]
    np.testing.assert_array_equal(result.layers["b4_soupx"], a.X)


def test_report_lists_excluded_samples(tmp_path, monkeypatch):
    report = load_module("report", "07_report.py")
    a = counts(np.ones((3, 5)))
    a.obs["Sample"] = ["included"] * 3
    a.uns["sample_inclusion"] = pd.DataFrame({
        "sample": ["included", "excluded"],
        "status": ["included", "excluded"],
        "reason": ["", "fewer than 3 nonzero-count cells"],
    })
    path = tmp_path / "report_input.h5ad"
    a.write_h5ad(path)
    monkeypatch.setattr(report, "_plot_qc_violin", lambda a: None)
    html = report._section_qc(ad.read_h5ad(path))
    assert "excluded" in html
    assert "fewer than 3 nonzero-count cells" in html


def test_single_cluster_diagnostics_and_umap_failure(tmp_path, monkeypatch):
    a = counts(np.random.default_rng(3).poisson(2, (10, 40)))
    a.obs["Sample"] = "A"

    def single_cluster(adata, key_added, **kwargs):
        adata.obs[key_added] = pd.Categorical(["0"] * adata.n_obs)

    def no_umap(*args, **kwargs):
        raise np.linalg.LinAlgError("embedding failed")

    monkeypatch.setattr(qc.sc.tl, "leiden", single_cluster)
    monkeypatch.setattr(qc.sc.tl, "umap", no_umap)
    result = qc.cluster_and_embed(a, [1.0])
    assert result.obs["leiden_1_0"].nunique() == 1
    assert "rank_genes_groups_1_0" not in result.uns
    assert result.uns["diagnostic_failures"]["umap"] == "embedding failed"
    qc.save_qc_plots(result, "A", str(tmp_path), [1.0])
