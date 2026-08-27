"""Continuity heuristics + Palantir/DPT/scVelo wiring (packages optional)."""

from __future__ import annotations

from pathlib import Path

import pytest

from agents.artifacts import FIGURE_CAPTIONS
from agents.code_schema import validate_script
from agents.intent import rule_intent
from agents.templates import cluster_annotate_script
from rag.synonyms import expand_query
from scagent.trajectory import (
    inspect_trajectory_hints,
    plot_gene_trends,
    should_plan_trajectory,
)


def test_pbmc_hints_not_candidate():
    hints = inspect_trajectory_hints(tissue="pbmc", n_obs=3000)
    assert hints["trajectory_candidate"] is False
    assert hints["has_velocity_layers"] is False


def test_embryo_and_velocity_layers_are_candidates():
    embryo = inspect_trajectory_hints(tissue="embryo", n_obs=800)
    assert embryo["trajectory_candidate"] is True
    vel = inspect_trajectory_hints(tissue="pbmc", layers=["spliced", "unspliced"], n_obs=800)
    assert vel["trajectory_candidate"] is True
    assert vel["has_velocity_layers"] is True


def test_user_query_overrides_discrete_tissue():
    hints = inspect_trajectory_hints(tissue="pbmc", query="用 Palantir 看命运", n_obs=800)
    assert hints["trajectory_candidate"] is True


def test_should_plan_trajectory_respects_config_and_tissue():
    pbmc = {"tissue": "pbmc", "n_cells": 2000}
    assert should_plan_trajectory(pbmc, {"intents": ["qc", "annotation"]}) is False
    assert should_plan_trajectory({"tissue": "hematopoiesis"}, {"intents": ["qc"]}) is True
    assert should_plan_trajectory(pbmc, {"intents": ["qc"]}, {"modules": {"trajectory": "off"}}) is False
    assert should_plan_trajectory(pbmc, {"intents": ["qc"]}, {"modules": {"trajectory": "force"}}) is True
    assert should_plan_trajectory({"trajectory_candidate": True}, {"intents": ["qc"]}) is True


def test_rule_intent_palantir_fate():
    intent = rule_intent("Palantir 命运")
    assert "trajectory" in intent["intents"]
    assert "trajectory" in rule_intent("scVelo RNA velocity")["intents"]


def test_downstream_template_calls_run_trajectory_phase():
    down = cluster_annotate_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"},
        {"nmads": 5},
        {"integrator": None},
    )
    assert "run_trajectory_phase" in down
    assert "sc.tl.dpt(" not in down
    r = validate_script(down, phase="downstream")
    assert r["ok"] is True, r["issues"]
    assert r["steps"].index("leiden") < r["steps"].index("trajectory")


def test_caption_prefers_pseudotime_over_umap():
    low = "umap_pseudotime.png"
    kind = "other"
    for key in sorted(FIGURE_CAPTIONS, key=len, reverse=True):
        if key in low:
            kind = key
            break
    assert kind == "pseudotime"
    assert "gene_trends" in FIGURE_CAPTIONS
    assert "paga" in FIGURE_CAPTIONS
    assert "velocity" in FIGURE_CAPTIONS


def test_synonym_expands_trajectory():
    q = expand_query("轨迹分析").lower()
    assert "palantir" in q
    assert "scvelo" in q
    assert "monocle3" in q


def test_dual_downstream_mentions_trajectory_verdict():
    from scagent.dual import render_dual_markdown

    md = render_dual_markdown(
        {
            "user_query": "命运",
            "report_lang": "zh",
            "code_downstream": "from scagent.trajectory import run_trajectory_phase\n",
            "execution_downstream": {"executed": True, "ok": True},
            "artifacts": {
                "metrics": {"trajectory_verdict": "run", "trajectory_methods": ["paga", "dpt"]},
            },
            "plan": {"route": ["trajectory"]},
        }
    )
    assert "轨迹" in md
    assert "paga" in md


def test_assess_continuity_skips_tiny_or_discrete():
    pytest.importorskip("anndata")
    pytest.importorskip("numpy")
    import anndata as ad
    import numpy as np
    import pandas as pd

    from scagent.trajectory import assess_continuity

    tiny = ad.AnnData(np.ones((10, 5)))
    out = assess_continuity(tiny)
    assert out["verdict"] == "skip"
    forced = assess_continuity(tiny, force=True)
    assert forced["verdict"] == "run"
    assert forced.get("confidence") == "low"

    X = np.random.default_rng(0).poisson(2, size=(60, 20)).astype(float)
    adata = ad.AnnData(X)
    adata.obs["leiden"] = pd.Categorical(["0"] * 20 + ["1"] * 20 + ["2"] * 20)
    adata.obsm["X_umap"] = np.column_stack(
        [np.repeat([0.0, 5.0, 10.0], 20), np.zeros(60)]
    )
    blob = assess_continuity(adata)
    assert "n_clusters" in blob
    assert blob["verdict"] in {"run", "skip"}


def test_run_trajectory_phase_force_on_synthetic(tmp_path: Path):
    pytest.importorskip("scanpy")
    pytest.importorskip("anndata")
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc

    from scagent.trajectory import run_trajectory_phase

    rng = np.random.default_rng(1)
    n = 80
    t = np.linspace(0, 1, n)
    X = rng.poisson(np.outer(1 + 3 * t, np.linspace(0.5, 2.0, 25))).astype(float)
    adata = ad.AnnData(X)
    adata.obs["leiden"] = pd.Categorical(pd.cut(t, 4, labels=["0", "1", "2", "3"]).astype(str))
    adata.obs_names = [f"c{i}" for i in range(n)]
    adata.var_names = [f"g{i}" for i in range(adata.n_vars)]
    adata.var["highly_variable"] = True
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=8)
    sc.pp.neighbors(adata, n_neighbors=8, n_pcs=8)
    sc.tl.umap(adata)
    payload = run_trajectory_phase(adata, mode="force", workspace=tmp_path)
    assert payload.get("verdict") == "run"
    assert "dpt" in (payload.get("methods") or []) or "paga" in (payload.get("methods") or [])
    assert "dpt_pseudotime" in adata.obs
    assert (tmp_path / "figures" / "gene_trends.png").is_file() or payload.get("gene_trends")


def test_plot_gene_trends_writes_png(tmp_path: Path):
    pytest.importorskip("anndata")
    pytest.importorskip("matplotlib")
    import anndata as ad
    import numpy as np

    adata = ad.AnnData(np.random.default_rng(2).random((40, 6)))
    adata.obs["dpt_pseudotime"] = np.linspace(0, 1, 40)
    adata.var_names = [f"GENE{i}" for i in range(6)]
    path = plot_gene_trends(adata, figdir=tmp_path)
    assert path
    assert Path(path).is_file()
