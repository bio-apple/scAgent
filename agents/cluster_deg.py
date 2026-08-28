"""Clustering & Differential Agent: Leiden, annotation evidence, cluster markers, condition DEG."""

from __future__ import annotations

from agents.annotation import build_annotation_plan
from scagent.best_practices_loader import practices_for_phase


def build_cluster_deg_plan(state: dict) -> dict:
    """Domain plan for clustering + DEG. Annotation stays here because DEG groupby needs cell types."""
    ann = build_annotation_plan(state)
    plan = state.get("plan") or {}
    return {
        "role": "cluster_deg",
        "agent": "Clustering & Differential Agent",
        "tasks": ["neighbors", "umap", "leiden", "annotation", "cluster_markers", "deg", "trajectory"],
        "annotation": ann,
        "best_practices": practices_for_phase(
            "downstream",
            route=list(plan.get("route") or []),
            query=state.get("user_query"),
        ),
        "needs_pseudobulk": bool(plan.get("needs_pseudobulk")),
        "force_pseudobulk_de": bool(plan.get("force_pseudobulk_de")),
        "n_replicates": plan.get("n_replicates"),
        "deg_engine": plan.get("deg_engine") or "auto",
        "marker_method": plan.get("marker_method") or "auto",
        "deg_cross_validate": plan.get("deg_cross_validate", "auto"),
        "condition_key": plan.get("condition_key"),
        "integrator": plan.get("integrator"),
        "code": ann.get("code") or "",
    }
