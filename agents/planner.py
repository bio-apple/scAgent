from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from rag.retriever import format_hits, retrieve
from scagent.config import analysis_params, load_config
from scagent.skills_loader import recommend_skills, skill_catalog_text

SCVI_CELL_CUTOFF = 100_000
SCVI_SAMPLE_CUTOFF = 8


def choose_integrator(meta: dict, requested: str | None = None) -> str | None:
    req = (requested or meta.get("integrator_requested") or "auto") or "auto"
    req = str(req).lower()
    if req in {"none", "off", "skip"}:
        return None
    if req in {"harmony", "scvi", "cca"}:
        return req
    n_samples = int(meta.get("n_samples") or 1)
    n_cells = int(meta.get("n_cells") or 0)
    if not meta.get("need_batch_correction") and n_samples <= 1:
        return None
    if n_cells >= SCVI_CELL_CUTOFF or n_samples >= SCVI_SAMPLE_CUTOFF:
        return "scvi"
    return "harmony"


def build_plan(state: dict) -> dict:
    meta = dict(state.get("metadata") or {})
    language = state.get("language") or "python"
    if state.get("batch_key"):
        meta["sample_key"] = state["batch_key"]
        meta["need_batch_correction"] = True
    if state.get("resolution") is not None:
        meta["resolution"] = state["resolution"]
    if state.get("markers_path"):
        meta["markers_path"] = state["markers_path"]

    r_degraded = language == "r"
    cfg = load_config()
    requested = state.get("integrator") or (cfg.get("modules") or {}).get("batch") or "auto"
    integrator = None if r_degraded else choose_integrator(meta, requested)
    if str(requested).lower() in {"none", "off", "skip"}:
        meta["skip_integration_reason"] = "user disabled batch module"
    skills = [] if r_degraded else recommend_skills({**meta, "integrator": integrator}, language=language)
    imputation = "none" if r_degraded else (state.get("imputation") or (cfg.get("modules") or {}).get("imputation") or "none")
    resolution = state.get("resolution")
    if resolution is None:
        resolution = analysis_params(cfg).get("leiden_resolution")
    rag = format_hits(
        retrieve(
            f"{meta.get('platform')} {meta.get('tissue')} scRNA-seq best practices integration Harmony scVI",
            collections=["papers", "best_practices", "methods"],
        )
    )
    if r_degraded:
        route = ["plan_only"]
    else:
        route = ["qc", "normalize", "hvg", "pca"]
        if integrator:
            route.append(integrator)
        if imputation and imputation != "none":
            route.append(f"impute_{imputation}")
        route += ["neighbors", "leiden", "umap", "annotate"]
        if "deg" in str(state.get("user_query") or "").lower() or "差异" in str(state.get("user_query") or ""):
            route.append("pseudobulk_deg")
        if "轨迹" in str(state.get("user_query") or "") or "paga" in str(state.get("user_query") or "").lower():
            route.append("trajectory")

    plan = {
        "objective": state.get("user_query") or "标准 scRNA-seq 分析",
        "species": meta.get("species"),
        "platform": meta.get("platform"),
        "n_samples": meta.get("n_samples"),
        "n_cells": meta.get("n_cells"),
        "need_batch_correction": bool(meta.get("need_batch_correction")),
        "integrator": integrator,
        "skip_integration_reason": meta.get("skip_integration_reason"),
        "imputation": imputation or "none",
        "resolution": resolution,
        "language": language,
        "r_degraded": r_degraded,
        "skills": skills,
        "route": route,
        "risks": [],
        "rag_excerpt": rag,
    }
    if r_degraded:
        plan["risks"].append(
            "现有 skills 为 Python/Scanpy。--language r 已降级为仅规划，不生成半成品 Seurat 代码。"
        )
    if meta.get("platform") == "parse":
        plan["risks"].append("Parse 平台 barcode 与 10x 不同，加载时不要套 Cell Ranger 默认假设。")
    if integrator:
        plan["risks"].append(
            f"多样本选用 {integrator}（n_cells={meta.get('n_cells')}, n_samples={meta.get('n_samples')}）。"
            "禁止把 UMAP 混匀当整合成功。"
        )
    if language != "r":
        llm = run_specialist(
            read_prompt("planner"),
            (
                f"用户任务: {state.get('user_query')}\n"
                f"metadata: {json.dumps(meta, ensure_ascii=False)}\n"
                f"integrator={integrator}\n"
                f"可选 skills:\n{skill_catalog_text()}\n"
                f"RAG:\n{rag}\n"
                "请输出分析路线。R/Seurat 不要生成代码。"
            ),
        )
        plan["narrative"] = llm or (
            f"物种={plan['species']}，平台={plan['platform']}，样本数={plan['n_samples']}。"
            f"整合={integrator or '无'}。路线: {' → '.join(route)}。skills: {', '.join(skills)}。"
        )
    else:
        plan["narrative"] = (
            "R/Seurat 路径未实现为可执行代码（仓库 SOP 为 Scanpy）。以下为规划与风险，请在 R 环境按 skills 意图手工实现，"
            "或改用 --language python。"
            f" 建议路线: QC(MAD) → 归一化 → HVG → PCA → 整合(如需) → Leiden → 注释双验证。"
        )
    return plan
