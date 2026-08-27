import pytest

from agents.templates import qc_preprocess_script


def _tiny_h5ad(path, n_cells=80, n_genes=40):
    import numpy as np
    from anndata import AnnData

    rng = np.random.default_rng(0)
    X = rng.poisson(1.5, size=(n_cells, n_genes)).astype(np.float32)
    var_names = [f"Gene{i}" for i in range(n_genes)]
    var_names[0] = "MT-ND1"
    var_names[1] = "MT-CO1"
    var_names[2] = "CD3D"
    obs_names = [f"c{i}" for i in range(n_cells)]
    adata = AnnData(X)
    adata.obs_names = obs_names
    adata.var_names = var_names
    adata.obs["sample"] = ["s1" if i < n_cells // 2 else "s2" for i in range(n_cells)]
    adata.write_h5ad(path)
    return path


def test_qc_script_on_synthetic_h5ad(tmp_path):
    anndata = pytest.importorskip("anndata")
    pytest.importorskip("scanpy")
    del anndata
    h5ad = _tiny_h5ad(tmp_path / "tiny.h5ad")
    code = qc_preprocess_script(
        {
            "data_path": str(h5ad),
            "species": "human",
            "tissue": "pbmc",
            "sample_key": "sample",
            "need_batch_correction": True,
        },
        {"nmads": 5, "extra_qc": ["hb"]},
    )
    compile(code, "<qc_syn>", "exec")
    from sandbox.executor import write_and_maybe_run

    result = write_and_maybe_run(
        code,
        workspace=tmp_path,
        execute=True,
        timeout=120,
        filename="qc_preprocess.py",
        extra_manifest={"data_path": str(h5ad)},
    )
    if not result["ok"] and result.get("missing_packages"):
        pytest.skip("analysis extras missing: " + ",".join(result["missing_packages"]))
    assert result["ok"], result["stderr"][-800:]
    assert (tmp_path / "adata_qc.h5ad").exists()
    assert "SCAGENT_METRICS:" in result["stdout"]
    assert any("violin" in p.lower() or "scatter" in p.lower() for p in result["figures"])
