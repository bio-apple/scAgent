import numpy as np
import pytest

from agents.qc_expert import build_qc_strategy
from agents.reviewer import audit_code, publication_review
from agents.templates import qc_preprocess_script
from scagent.demo import write_tiny_h5ad
from scagent.doublets import (
    _consensus,
    _simulate_doublet_scores,
    detect_doublets,
    expected_doublet_rate,
    resolve_doublet_methods,
)


def test_resolve_doublet_methods_auto():
    assert resolve_doublet_methods("auto", tissue="pbmc", n_samples=1) == ["scrublet"]
    assert resolve_doublet_methods("auto", tissue="pbmc", n_samples=3) == ["scrublet", "scdblfinder"]
    assert resolve_doublet_methods("auto", tissue="tumor", n_samples=1) == ["scrublet", "scdblfinder"]
    assert resolve_doublet_methods("auto", tissue="brain", n_samples=1) == ["scrublet", "scdblfinder"]
    assert resolve_doublet_methods("scrublet", tissue="tumor", n_samples=8) == ["scrublet"]
    assert resolve_doublet_methods("both", tissue="pbmc", n_samples=1) == ["scrublet", "scdblfinder"]


def test_expected_rate_10x_rule():
    assert expected_doublet_rate(1000) == pytest.approx(0.008)
    assert expected_doublet_rate(20_000) == pytest.approx(0.10)


def test_consensus_intersection_and_discordant():
    a = np.array([True, True, False, False])
    b = np.array([True, False, True, False])
    both, disc, agree = _consensus(a, b)
    assert both.tolist() == [True, False, False, False]
    assert disc.tolist() == [False, True, True, False]
    assert agree == pytest.approx(0.5)


def test_count_simulation_on_tiny_h5ad(tmp_path):
    pytest.importorskip("sklearn")
    pytest.importorskip("anndata")
    path = write_tiny_h5ad(tmp_path / "tiny.h5ad")
    import anndata as ad

    adata = ad.read_h5ad(path)
    score, pred = _simulate_doublet_scores(adata)
    assert len(score) == adata.n_obs
    assert pred.dtype == bool
    assert 0 < int(pred.sum()) < adata.n_obs
    out = detect_doublets(adata, methods=["sim"], sample_key="sample", tissue="pbmc", remove=False)
    info = out.uns["doublets"]
    assert info["status"] in {"ok", "partial"}
    assert "sim" in info["methods"] or "scdblfinder" in info["methods"]
    assert "predicted_doublet" in out.obs
    assert "doublet_discordant" in out.obs


def test_qc_template_crosscheck_for_tumor():
    code = qc_preprocess_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "tumor"},
        {"nmads": 6, "doublet_methods": "auto"},
    )
    compile(code, "<dbl>", "exec")
    assert "detect_doublets" in code
    assert "methods='both'" in code or 'methods="both"' in code
    assert "scDblFinder" in code
    r = audit_code(code, {"tissue": "tumor", "n_samples": 1}, phase="qc")
    assert r["passed"] is True
    assert not any(i.get("id") == "qc.doublet_cross" for i in r.get("issue_records") or [])


def test_qc_strategy_multi_sample_requests_crosscheck():
    s = build_qc_strategy({"metadata": {"tissue": "pbmc", "n_samples": 4, "need_batch_correction": True}})
    assert s["doublet_methods_resolved"] == ["scrublet", "scdblfinder"]
    assert "交叉验证" in (s.get("protocol") or "")


def test_publication_mentions_crosscheck_when_both_engines():
    card = publication_review(
        {
            "code_qc": "from scagent.doublets import detect_doublets\npredicted_doublet\nmethods='both' scDblFinder",
            "review_qc": {"passed": True},
            "artifacts": {
                "metrics": {
                    "doublet_status": "ok",
                    "doublet_rate": 0.03,
                    "doublet_agreement": 0.91,
                    "doublet_methods": ["scrublet", "sim"],
                }
            },
            "metadata": {"n_samples": 1},
        }
    )
    by_key = {i["key"]: i for i in card["items"]}
    assert by_key["doublet_detection"]["status"] == "pass"
    assert "交叉验证" in by_key["doublet_detection"]["detail"]
