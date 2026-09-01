import pytest

from agents.markers import choose_celltypist_model
from agents.qc_expert import build_qc_strategy
from agents.reviewer import audit_code
from agents.templates import cluster_annotate_script, qc_preprocess_script
from scagent.analysis import knn_ilisi, needs_condition_de, pca_batch_r2
from scagent.preprocess import choose_ambient


def test_celltypist_model_by_tissue():
    assert choose_celltypist_model("pbmc") == "Immune_All_Low.pkl"
    assert choose_celltypist_model("lung") == "Human_Lung_Atlas.pkl"
    assert choose_celltypist_model("brain") == "Developing_Human_Brain.pkl"
    assert choose_celltypist_model("heart") == "Adult_Human_Heart.pkl"
    assert choose_celltypist_model("liver") == "Adult_Human_Liver.pkl"
    assert choose_celltypist_model("default") is None


def test_qc_writes_predicted_doublet_and_cell_cycle():
    code = qc_preprocess_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"},
        {"nmads": 5, "remove_doublets": True, "regress_cell_cycle": "auto"},
    )
    compile(code, "<qc_adv>", "exec")
    assert "detect_doublets" in code
    assert "doublet_call" in code or "DOUBLET_FILTER" in code
    assert "REMOVE_DOUBLETS = True" in code
    assert "cell_cycle_score" in code
    assert "scrublet skipped" not in code
    r = audit_code(code, {"tissue": "pbmc"}, phase="qc")
    assert r["passed"] is True


def test_brain_ambient_is_real_correction():
    assert choose_ambient("brain", "auto") == "soupx"
    assert choose_ambient("pbmc", "auto") == "none"
    s = build_qc_strategy({"metadata": {"tissue": "brain"}})
    assert s["ambient"] == "soupx"
    code = qc_preprocess_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "brain"},
        s,
    )
    compile(code, "<amb>", "exec")
    assert "remove_ambient" in code
    assert "consider SoupX/DecontX" not in code


def test_condition_de_requires_pseudobulk_impl():
    assert needs_condition_de("比较对照组 vs 处理组的差异表达") is True
    assert needs_condition_de("对 PBMC 做标准注释") is False
    assert needs_condition_de("degradation analysis") is False
    meta = {
        "data_path": "x.h5ad",
        "species": "human",
        "need_batch_correction": True,
        "sample_key": "batch",
        "tissue": "pbmc",
        "needs_pseudobulk": True,
        "condition_key": "condition",
        "n_replicates": 2,
        "force_pseudobulk_de": True,
    }
    code = cluster_annotate_script(
        meta,
        {},
        {
            "integrator": "harmony",
            "needs_pseudobulk": True,
            "condition_key": "condition",
            "force_pseudobulk_de": True,
        },
    )
    compile(code, "<pb>", "exec")
    assert "pseudobulk_de" in code
    assert "engine=" in code
    r = audit_code(code, meta, phase="downstream")
    assert r["has_pseudobulk_impl"] is True
    assert r["passed"] is True
    bad = code.replace("pseudobulk_de(", "wilcoxon_only(")
    r2 = audit_code(bad, meta, phase="downstream")
    assert r2["passed"] is False
    assert any("pseudobulk" in x for x in r2["issues"])


def test_liver_does_not_use_immune_all():
    meta = {"data_path": "x.h5ad", "tissue": "liver", "species": "human"}
    code = cluster_annotate_script(meta, {}, {"integrator": None})
    compile(code, "<liver>", "exec")
    assert "Immune_All_Low.pkl" not in code
    assert "Adult_Human_Liver.pkl" in code
    assert "ref2_label" in code
    r = audit_code(code, meta, phase="downstream")
    assert r["passed"] is True
    hijack = code.replace("Adult_Human_Liver.pkl", "Immune_All_Low.pkl")
    r2 = audit_code(hijack, meta, phase="downstream")
    assert r2["passed"] is False


def test_integration_metrics_helpers():
    np = pytest.importorskip("numpy")
    pytest.importorskip("sklearn")

    rng = np.random.default_rng(0)
    x = rng.normal(size=(40, 5))
    batch = ["a"] * 20 + ["b"] * 20
    r2 = pca_batch_r2(x, batch)
    assert 0.0 <= r2 <= 1.0
    mixed = knn_ilisi(x, batch, k=5)
    assert 0.0 <= mixed <= 1.0
    down = cluster_annotate_script(
        {"data_path": "x.h5ad", "need_batch_correction": True, "sample_key": "batch", "tissue": "pbmc"},
        {},
        {"integrator": "harmony"},
    )
    assert "integration_quality" in down
    assert "ilisi" in down
    assert "integration_plots" in down
    assert "integ_plots" in down
