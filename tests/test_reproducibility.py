from pathlib import Path

import pytest

from sandbox.executor import write_manifest
from scagent.reproducibility import (
    build_step_provenance,
    compute_environment_fingerprint,
    seed_propagation_record,
    summarize_adata,
)


def test_environment_fingerprint_has_hash():
    fp = compute_environment_fingerprint()
    assert len(fp["hash"]) == 64
    assert "package_versions" in fp["sources"]
    assert fp["packages"]["scagent"]


def test_seed_propagation_lists_steps():
    rec = seed_propagation_record(42)
    assert rec["master_seed"] == 42
    assert rec["steps"]["leiden"] == 42
    assert rec["steps"]["umap"] == 42
    assert rec["steps"]["hvg"] == 42


def test_write_manifest_includes_environment_and_seed(tmp_path: Path):
    write_manifest(tmp_path, {"phase": "qc", "data_path": "/tmp/x.h5ad"})
    import json

    payload = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert payload["scagent_version"]
    assert payload["environment"]["hash"]
    assert payload["seed_propagation"]["master_seed"] is not None
    assert "steps" in payload["seed_propagation"]


def test_summarize_adata_minimal():
    ad = pytest.importorskip("anndata")
    import numpy as np

    a = ad.AnnData(np.zeros((5, 3)))
    a.obs["sample"] = ["a"] * 5
    a.var["mt"] = [False, True, False]
    s = summarize_adata(a, step="test")
    assert s["n_obs"] == 5
    assert s["n_vars"] == 3
    assert "sample" in s["obs_columns"]
    assert "mt" in s["var_columns"]


def test_build_step_provenance_qc(tmp_path: Path, monkeypatch):
    ad = pytest.importorskip("anndata")
    import numpy as np

    a = ad.AnnData(np.zeros((10, 4)))
    a.write(tmp_path / "adata_qc.h5ad")
    steps = build_step_provenance(tmp_path, phase="qc", data_path=str(tmp_path / "adata_qc.h5ad"))
    assert len(steps) == 1
    assert steps[0]["step"] == "qc"
    assert steps[0]["output"]["n_obs"] == 10


def test_dask_params_default_disabled():
    from scagent.config import dask_params, gpu_params, load_config

    load_config(reload=True)
    assert dask_params()["enabled"] is False
    assert gpu_params()["enabled"] is False


def test_analysis_random_state_propagation():
    pytest.importorskip("scanpy")
    ad = pytest.importorskip("anndata")
    import numpy as np
    import scanpy as sc

    from scagent.analysis import leiden, neighbors, pca, umap

    rng = np.random.default_rng(0)
    x = rng.poisson(1, size=(40, 20)).astype(float)
    a = ad.AnnData(x)
    a.obs["leiden"] = ["0"] * 20 + ["1"] * 20
    a.var["highly_variable"] = True

    sc.pp.log1p(a)
    sc.pp.scale(a)
    sc.tl.pca(a, n_comps=5)
    pca(a, n_pcs=5, random_state=7)
    neighbors(a, n_neighbors=5, n_pcs=5, random_state=7)
    umap(a, random_state=7)
    leiden(a, resolution=0.5, random_state=7)
    assert "X_umap" in a.obsm


def test_report_includes_manifest_when_present(tmp_path, monkeypatch):
    from agents.writer import render_report

    ws = tmp_path / "workspace"
    ws.mkdir()
    import json

    (ws / "run_manifest.json").write_text(
        json.dumps(
            {
                "seed": 0,
                "environment": {"hash": "abc123", "sources": ["pip_freeze"]},
                "seed_propagation": {"master_seed": 0, "steps": {"leiden": 0, "umap": 0}},
                "step_provenance": [
                    {
                        "step": "qc",
                        "input": {"n_obs": 1000},
                        "output": {"n_obs": 900, "obs_columns": ["sample", "pct_counts_mt"]},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr("scagent.config.resolve_path", lambda cfg, key: ws if key == "workspace" else tmp_path)
    md = render_report(
        {
            "user_query": "t",
            "report_lang": "zh",
            "metadata": {},
            "plan": {},
            "artifacts": {},
            "review_publication": {"items": [], "score": 0, "max_score": 100, "verdict": "PASS", "passed": True},
        }
    )
    assert "abc123" in md
    assert "步骤 I/O" in md
