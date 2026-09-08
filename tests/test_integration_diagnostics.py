from types import SimpleNamespace

import anndata as ad
import numpy as np
import pytest

from test_sample_eligibility import counts, load_module


def test_failed_pacmap_cannot_destroy_core_output(tmp_path, monkeypatch):
    diagnostics = load_module("integration_diagnostics", "03_integration_diagnostics.py")
    a = counts(np.random.default_rng(2).poisson(2, (20, 30)))
    a.obsm["X_pca_harmony"] = np.ones((20, 5))
    source = tmp_path / "integrated.h5ad"
    a.write_h5ad(source)
    before = source.read_bytes()

    class BrokenPaCMAP:
        def __init__(self, **kwargs):
            pass

        def fit_transform(self, x):
            raise RuntimeError("PaCMAP failure")

    monkeypatch.setattr(diagnostics.pacmap, "PaCMAP", BrokenPaCMAP)
    with pytest.raises(RuntimeError, match="PaCMAP failure"):
        diagnostics.main(SimpleNamespace(
            input=str(source), output=str(tmp_path / "diagnostics/out.h5ad"),
            leiden_resolutions="1.0",
        ))
    assert source.read_bytes() == before
    assert ad.read_h5ad(source).n_obs == 20


def test_diagnostics_cannot_overwrite_integration(tmp_path):
    diagnostics = load_module("integration_diagnostics", "03_integration_diagnostics.py")
    with pytest.raises(ValueError, match="must not overwrite"):
        diagnostics.main(SimpleNamespace(input=str(tmp_path / "same.h5ad"),
                                         output=str(tmp_path / "same.h5ad")))
