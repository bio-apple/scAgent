"""P0/P1 audit follow-ups: ambient honesty, tissue unannotated, per-sample QC, resolution joint."""

from __future__ import annotations

import numpy as np
import pytest

from agents.templates import cluster_annotate_script, qc_preprocess_script
from scagent.analysis import choose_leiden_resolution
from scagent.preprocess import choose_ambient, filter_dynamic


def test_choose_ambient_auto_is_none_until_backend_wired():
    assert choose_ambient("brain", "auto") == "none"
    assert choose_ambient("tumor", "auto") == "none"
    assert choose_ambient("pbmc", "auto") == "none"
    assert choose_ambient("brain", "soupx") == "soupx"
    assert choose_ambient("brain", "none") == "none"


def test_qc_template_has_per_sample_and_normalization():
    code = qc_preprocess_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc", "sample_key": "sample"},
        {"nmads": 5, "normalization": "log1p", "per_sample_qc": True},
    )
    compile(code, "<qc_p1>", "exec")
    assert "PER_SAMPLE_QC" in code
    assert "normalize_expression" in code
    assert "log1p" in code


def test_cluster_template_unannotated_when_tissue_unknown():
    code = cluster_annotate_script(
        {"data_path": "x.h5ad", "tissue": "unknown", "species": "human"},
        {},
        {"integrator": None},
    )
    compile(code, "<ann_skip>", "exec")
    assert "ANNOTATION_SKIPPED" in code
    assert "unannotated" in code
    assert "SCAGENT_ANNOTATION_STATUS" in code


def test_cluster_template_joint_resolution():
    code = cluster_annotate_script(
        {"data_path": "x.h5ad", "tissue": "pbmc", "species": "human"},
        {},
        {"integrator": None},
    )
    compile(code, "<res_joint>", "exec")
    assert "choose_leiden_resolution" in code
    assert "marker_interpretability_score" in code


def test_choose_leiden_resolution_penalizes_two_clusters():
    sil = {0.2: 0.55, 0.4: 0.50, 0.6: 0.48, 0.8: 0.45, 1.0: 0.42}
    ncl = {0.2: 2, 0.4: 4, 0.6: 6, 0.8: 9, 1.0: 11}
    chosen, detail = choose_leiden_resolution(sil, ncl, marker_scores=None, default=0.6)
    assert chosen != 0.2
    assert chosen in {0.4, 0.6, 0.8, 1.0}
    assert detail["reason"] == "joint_silhouette_marker_size"


def test_filter_dynamic_per_sample():
    pytest.importorskip("anndata")
    from anndata import AnnData
    from scipy import sparse

    rng = np.random.default_rng(0)
    n = 80
    X = sparse.csr_matrix(rng.poisson(5, size=(n, 30)).astype(float))
    adata = AnnData(X)
    adata.obs["sample"] = ["A"] * 40 + ["B"] * 40
    # Sample B has inflated mito
    adata.obs["pct_counts_mt"] = np.concatenate([rng.normal(2, 0.5, 40), rng.normal(20, 1.0, 40)])
    adata.obs["log1p_total_counts"] = np.log1p(rng.normal(1000, 50, n))
    adata.obs["log1p_n_genes_by_counts"] = np.log1p(rng.normal(500, 20, n))
    adata.obs["n_genes_by_counts"] = rng.integers(200, 800, n)
    adata.obs["total_counts"] = rng.integers(500, 2000, n)
    out = filter_dynamic(adata, method="mad", nmads=3, sample_key="sample", per_sample=True)
    assert out.n_obs < n
    assert "outlier" in adata.obs
