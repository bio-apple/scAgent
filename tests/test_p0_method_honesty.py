"""P0 method honesty: ambient no silent heuristic; confirmatory PB refuses t-test; R QC continues."""

from __future__ import annotations

import numpy as np
import pytest

from agents.reviewer import audit_code, publication_review
from agents.templates import cluster_annotate_script, qc_preprocess_script
from scagent.preprocess import remove_ambient


def _toy_counts(n_obs: int = 40, n_vars: int = 20):
    pytest.importorskip("anndata")
    from anndata import AnnData
    from scipy import sparse

    rng = np.random.default_rng(0)
    X = rng.poisson(2.0, size=(n_obs, n_vars)).astype(np.float32)
    adata = AnnData(sparse.csr_matrix(X))
    adata.obs_names = [f"c{i}" for i in range(n_obs)]
    adata.var_names = [f"G{i}" for i in range(n_vars)]
    return adata


def test_ambient_does_not_mutate_without_heuristic_flag():
    adata = _toy_counts()
    before = adata.X.copy()
    out = remove_ambient(adata, method="soupx", allow_heuristic=False)
    info = out.uns.get("ambient") or {}
    assert info.get("applied") is False
    assert "unavailable" in str(info.get("method"))
    assert (out.X != before).nnz == 0 if hasattr(out.X, "nnz") else np.allclose(out.X.toarray(), before.toarray())


def test_ambient_heuristic_only_when_allowed():
    adata = _toy_counts()
    before = adata.X.copy()
    out = remove_ambient(adata, method="soupx", allow_heuristic=True)
    info = out.uns.get("ambient") or {}
    assert info.get("applied") is True
    assert "heuristic" in str(info.get("method"))
    assert (out.X != before).nnz > 0


def test_qc_template_continues_after_r_qc():
    code = qc_preprocess_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "brain"},
        {"nmads": 5, "ambient": "soupx"},
    )
    compile(code, "<qc_p0>", "exec")
    assert "SCAGENT_R_QC_OK" in code
    assert "_R_QC" in code
    assert "raise SystemExit(0)" not in code
    assert "detect_doublets" in code
    assert "allow_heuristic=False" in code
    r = audit_code(code, {"tissue": "brain"}, phase="qc")
    assert r["passed"] is True


def test_confirmatory_pseudobulk_refuses_ttest_fallback(tmp_path, monkeypatch):
    from scagent.analysis import pseudobulk_de
    from tests.test_deg import _condition_adata

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("scagent.analysis._deg_via_rpy2", lambda *a, **k: None)
    monkeypatch.setattr("scagent.analysis._deg_via_rscript", lambda *a, **k: None)
    adata = _condition_adata()
    with pytest.raises(RuntimeError, match="confirmatory"):
        pseudobulk_de(
            adata,
            sample_key="sample",
            condition_key="condition",
            groupby="cell_type",
            min_cells=5,
            engine="auto",
            cross_validate=False,
            confirmatory=True,
        )


def test_forced_template_passes_confirmatory():
    code = cluster_annotate_script(
        {
            "data_path": "x.h5ad",
            "tissue": "pbmc",
            "condition_key": "condition",
            "force_pseudobulk_de": True,
            "n_replicates": 2,
            "sample_key": "sample",
        },
        {},
        {
            "needs_pseudobulk": True,
            "force_pseudobulk_de": True,
            "condition_key": "condition",
            "integrator": None,
        },
    )
    compile(code, "<pb_conf>", "exec")
    assert "confirmatory=True" in code


def test_publication_fails_ttest_when_force_pseudobulk():
    card = publication_review(
        {
            "metadata": {
                "tissue": "pbmc",
                "condition_key": "condition",
                "n_replicates": 3,
                "force_pseudobulk_de": True,
            },
            "plan": {"force_pseudobulk_de": True, "needs_pseudobulk": True},
            "code_downstream": "pseudobulk_de(adata, confirmatory=True)\n",
            "review_qc": {"passed": True},
            "review_downstream": {"passed": True, "has_dual": True},
            "execute_code": True,
            "artifacts": {"metrics": {"deg_engine": "ttest_bh", "annotation_dual_validation": True}},
        }
    )
    deg = next(i for i in card["items"] if i["key"] == "deg")
    assert deg["status"] == "fail"
