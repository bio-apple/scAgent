"""Clustering & Differential Agent: Leiden, annotation evidence, cluster markers, condition DEG."""

from __future__ import annotations

from agents.annotation import build_annotation_plan


def build_cluster_deg_plan(state: dict) -> dict:
    """Domain plan for clustering + DEG. Annotation stays here because DEG groupby needs cell types."""
    ann = build_annotation_plan(state)
    plan = state.get("plan") or {}
    return {
        "role": "cluster_deg",
        "agent": "Clustering & Differential Agent",
        "tasks": ["neighbors", "umap", "leiden", "annotation", "cluster_markers", "deg", "trajectory"],
        "annotation": ann,
        "needs_pseudobulk": bool(plan.get("needs_pseudobulk")),
        "deg_engine": plan.get("deg_engine") or "auto",
        "marker_method": plan.get("marker_method") or "auto",
        "deg_cross_validate": plan.get("deg_cross_validate", "auto"),
        "condition_key": plan.get("condition_key"),
        "integrator": plan.get("integrator"),
        "code": ann.get("code") or "",
    }
