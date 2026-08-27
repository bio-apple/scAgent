import pytest

from scagent.inspect_data import detect_expression_layer, inspect_data


def _counts_adata(n_cells=40, n_genes=25, peak=80):
    import numpy as np
    from anndata import AnnData

    rng = np.random.default_rng(0)
    X = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32)
    X[0, 0] = peak
    adata = AnnData(X)
    adata.obs_names = [f"c{i}" for i in range(n_cells)]
    names = [f"G{i}" for i in range(n_genes)]
    names[0] = "MT-ND1"
    adata.var_names = names
    return adata


def test_detect_counts_log1p_scaled():
    pytest.importorskip("anndata")
    pytest.importorskip("numpy")
    adata = _counts_adata()
    info = detect_expression_layer(adata)
    assert info["layer"] == "counts"
    assert info["x_max"] >= 20
    assert info["uns_log1p"] is False

    adata.uns["log1p"] = {"base": None}
    assert detect_expression_layer(adata)["layer"] == "log1p"

    adata2 = _counts_adata()
    adata2.X = adata2.X.astype("float32")
    adata2.X = adata2.X - adata2.X.mean(axis=0)
    info_s = detect_expression_layer(adata2)
    assert info_s["layer"] == "scaled"
    assert info_s["x_min"] < 0


def test_inspect_h5ad_records_layer(tmp_path):
    pytest.importorskip("anndata")
    import numpy as np

    p = tmp_path / "counts.h5ad"
    _counts_adata().write_h5ad(p)
    meta = inspect_data(str(p))
    assert meta["expression_layer"] == "counts"
    assert meta["x_max"] is not None
    assert meta["sparsity"] is not None

    logged = _counts_adata()
    x = np.asarray(logged.X, dtype=np.float64)
    logged.X = np.log1p(x / np.maximum(x.sum(axis=1, keepdims=True), 1) * 1e4)
    logged.uns["log1p"] = {"base": None}
    lp = tmp_path / "log.h5ad"
    logged.write_h5ad(lp)
    meta2 = inspect_data(str(lp))
    assert meta2["expression_layer"] == "log1p"
    assert meta2["uns_log1p"] is True
    assert any("log1p" in n for n in meta2["notes"])


def test_normalize_log1p_is_idempotent():
    pytest.importorskip("scanpy")
    import numpy as np

    from scagent.preprocess import normalize_log1p

    adata = _counts_adata()
    normalize_log1p(adata)
    assert "log1p" in adata.uns
    x1 = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    normalize_log1p(adata)
    x2 = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    assert np.allclose(x1, x2)
    assert detect_expression_layer(adata)["layer"] == "log1p"


def test_normalize_log1p_rejects_scaled():
    pytest.importorskip("scanpy")
    import scanpy as sc

    from scagent.preprocess import normalize_log1p

    adata = _counts_adata()
    sc.pp.scale(adata, max_value=10)
    with pytest.raises(ValueError, match="scaled"):
        normalize_log1p(adata)


def test_pca_skips_second_scale():
    pytest.importorskip("scanpy")
    import numpy as np
    import scanpy as sc

    from scagent.analysis import pca
    from scagent.preprocess import normalize_log1p

    adata = _counts_adata(n_cells=50, n_genes=40, peak=120)
    normalize_log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=20, subset=False)
    pca(adata, n_pcs=8)
    x1 = np.asarray(adata.X)
    pca_x = np.array(adata.obsm["X_pca"])
    pca(adata, n_pcs=8)
    assert np.allclose(x1, np.asarray(adata.X))
    assert np.allclose(pca_x, np.asarray(adata.obsm["X_pca"]))


def test_qc_template_guards_double_normalize():
    from agents.templates import cluster_annotate_script, qc_preprocess_script

    code = qc_preprocess_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"},
        {"nmads": 5},
    )
    compile(code, "<qc_layer>", "exec")
    assert "detect_expression_layer" in code
    assert "skip normalize_total/log1p" in code
    down = cluster_annotate_script(
        {"data_path": "x.h5ad", "tissue": "pbmc"},
        {},
        {"integrator": None},
    )
    compile(down, "<down_layer>", "exec")
    assert "skip scale" in down
    assert "n_comps=" in down
