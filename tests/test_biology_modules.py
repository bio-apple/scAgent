from agents.markers import load_marker_catalog
from agents.planner import choose_integrator
from agents.qc_expert import build_qc_strategy
from agents.reviewer import audit_code
from agents.templates import cluster_annotate_script, qc_preprocess_script


def test_hierarchical_catalog():
    cat = load_marker_catalog(tissue="pbmc")
    tex = next(t for t in cat["cell_types"] if t["name"] == "CD8 Tex")
    assert tex["lineage"] == ["Immune", "T cell", "CD8 T", "CD8 Tex"]
    assert len(tex["positive"]) >= 2 and tex["negative"]


def test_percentile_and_no_hardcoded_mito():
    code = qc_preprocess_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "heart"},
        {
            "method": "percentile",
            "nmads": 6,
            "percentile": {"n_genes_low": 2, "n_genes_high": 98, "pct_mt_high": 98},
            "hard": {"pct_mt": None},
        },
    )
    compile(code, "<pct>", "exec")
    assert "np.percentile" in code
    assert "QC_METHOD" in code
    assert "HARD_PCT_MT = None" in code
    assert "pct_counts_mt < 5" not in code
    r = audit_code(code, {}, phase="qc")
    assert r["passed"] is True


def test_optional_cca_and_magic():
    meta = {"data_path": "x.h5ad", "need_batch_correction": True, "sample_key": "batch", "tissue": "pbmc"}
    down = cluster_annotate_script(meta, {}, {"integrator": "cca"})
    compile(down, "<cca>", "exec")
    assert "scanorama" in down.lower()
    assert "cell_type_l1" in down
    assert "annotation_conflict" in down
    r = audit_code(down, meta, phase="downstream")
    assert r["passed"] is True
    qc = qc_preprocess_script(meta, {"nmads": 5, "imputation": "magic"})
    compile(qc, "<magic>", "exec")
    assert "magic" in qc.lower()
    alra = qc_preprocess_script(meta, {"nmads": 5, "imputation": "alra"})
    compile(alra, "<alra>", "exec")
    assert "randomized_svd" in alra


def test_choose_integrator_optional_none():
    assert choose_integrator({"n_samples": 4, "need_batch_correction": True}, "none") is None
    assert choose_integrator({"n_samples": 1, "need_batch_correction": False}, "cca") == "cca"
    assert choose_integrator({"n_samples": 3, "need_batch_correction": True}, "auto") == "harmony"


def test_qc_strategy_reads_config_method():
    s = build_qc_strategy({"metadata": {"tissue": "kidney"}, "qc_method": "hybrid", "imputation": "alra"})
    assert s["method"] == "hybrid"
    assert s["imputation"] == "alra"
    assert s["hard"].get("pct_mt") in (None, "null") or s["hard"].get("pct_mt") is None
