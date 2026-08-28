from agents.planner import build_plan
from scagent.best_practices_loader import (
    list_practices,
    load_practice_text,
    practices_catalog_text,
    practices_for_phase,
    practices_for_route,
)


def test_reference_sops_present():
    names = {p.name for p in list_practices()}
    expected = {
        "qc",
        "doublet-detection",
        "normalization",
        "feature-selection",
        "dimensionality-reduction",
        "clustering",
        "cell-annotation",
        "marker-genes",
        "integration",
        "pseudobulk-de",
        "pathway-enrichment",
        "trajectory",
    }
    assert expected <= names
    assert "README" not in names


def test_practices_for_qc_phase():
    names = practices_for_phase("qc")
    assert names[:2] == ["qc", "doublet-detection"]
    assert "normalization" in names
    text = load_practice_text("qc")
    assert "MAD" in text or "mad" in text.lower()


def test_practices_for_deg_and_trajectory_query():
    names = practices_for_route(
        ["qc", "leiden", "annotate", "deg", "trajectory"],
        ["deg", "trajectory"],
        "肿瘤 vs 对照 pseudobulk 和拟时序",
    )
    assert "pseudobulk-de" in names
    assert "trajectory" in names
    assert "clustering" in names or "cell-annotation" in names


def test_catalog_lists_all_sops():
    text = practices_catalog_text()
    assert "qc" in text
    assert "pseudobulk-de" in text
    assert "knowledge/best_practices" in text


def test_planner_attaches_best_practices():
    plan = build_plan(
        {
            "user_query": "标准 PBMC QC 聚类注释",
            "metadata": {"tissue": "pbmc", "n_samples": 1, "species": "human", "n_cells": 100},
            "language": "python",
        }
    )
    bp = plan.get("best_practices") or []
    assert "qc" in bp
    assert "doublet-detection" in bp
    rag = plan.get("rag_excerpt") or ""
    assert rag
    assert "qc.md" in rag.lower() or "best_practices" in rag.lower() or "MAD" in rag
