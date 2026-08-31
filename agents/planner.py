from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from agents.dependencies import resolve_route, serialize_dag
from agents.intent import parse_intent
from agents.markers import choose_celltypist_model
from agents.roles import assign_roles
from rag.retriever import format_hits, retrieve_fused
from scagent.config import analysis_params, load_config
from scagent.deg_methods import force_pseudobulk_de, parse_deg_preference, resolve_forced_deg_engine
from scagent.preprocess import choose_ambient
from scagent.best_practices_loader import practices_catalog_text, practices_for_route
from scagent.skills_loader import recommend_skills, skill_catalog_text
from scagent.tool_router import analysis_language, build_tool_route

SCVI_CELL_CUTOFF = 100_000
SCVI_SAMPLE_CUTOFF = 8


def choose_integrator(meta: dict, requested: str | None = None) -> str | None:
    req = (requested or meta.get("integrator_requested") or "auto") or "auto"
    req = str(req).lower()
    if req in {"scanorama"}:
        req = "cca"
    if req in {"none", "off", "skip"}:
        return None
    if req in {"harmony", "scvi", "cca", "bbknn"}:
        return req
    n_samples = int(meta.get("n_samples") or 1)
    n_cells = int(meta.get("n_cells") or 0)
    if not meta.get("need_batch_correction") and n_samples <= 1:
        return None
    if meta.get("batch_condition_confounded"):
        return None
    if n_cells >= SCVI_CELL_CUTOFF or n_samples >= SCVI_SAMPLE_CUTOFF:
        return "scvi"
    return "harmony"


def explain_integrator(meta: dict, requested: str | None, chosen: str | None) -> str:
    """Human-readable reason for auto vs skip vs user override. Goes into the report."""
    req = str(requested or meta.get("integrator_requested") or "auto").lower()
    if req in {"scanorama"}:
        req = "cca"
    key = meta.get("sample_key") or "sample"
    n = int(meta.get("n_samples") or 1)
    n_cells = meta.get("n_cells")
    if req in {"none", "off", "skip"}:
        return "用户关闭批次模块（--integrator none / modules.batch=none）"
    if meta.get("batch_condition_confounded") and req == "auto":
        return f"obs[{key!r}] 与条件 1:1 共线，auto 跳过整合，以免把处理效应当批次抹掉"
    if chosen is None:
        if n <= 1 and not meta.get("need_batch_correction"):
            return f"未检测到多样本/批次列（sample_key={key}），不做整合"
        return str(meta.get("skip_integration_reason") or "未整合")
    if req in {"harmony", "scvi", "cca", "bbknn"}:
        label = {"cca": "Scanorama (cca)"}.get(chosen, chosen)
        return f"用户指定 {label}；批次列={key}，n_samples={n}，n_cells={n_cells}"
    if chosen == "scvi":
        return (
            f"检测到批次列 {key}（n_samples={n}，n_cells={n_cells}）："
            f"细胞数≥{SCVI_CELL_CUTOFF} 或样本数≥{SCVI_SAMPLE_CUTOFF}，auto 选 scVI"
        )
    if chosen == "harmony":
        return (
            f"检测到批次列 {key}（n_samples={n}，n_cells={n_cells}）："
            "auto 选 Harmony（Luecken 2022 简单–中等批次默认；--integrator bbknn|cca 可改）"
        )
    return f"integrator={chosen}；批次列={key}，n_samples={n}"


def build_plan(state: dict) -> dict:
    meta = dict(state.get("metadata") or {})
    cfg = load_config()
    language = state.get("language") or analysis_language(cfg)
    if state.get("batch_key"):
        meta["sample_key"] = state["batch_key"]
        meta["need_batch_correction"] = True
    if state.get("resolution") is not None:
        meta["resolution"] = state["resolution"]
    if state.get("markers_path"):
        meta["markers_path"] = state["markers_path"]

    r_degraded = language == "r"
    requested = state.get("integrator") or (cfg.get("modules") or {}).get("batch") or "auto"
    req_l = str(requested).lower()
    if meta.get("batch_condition_confounded") and req_l not in {"harmony", "scvi", "cca", "bbknn", "scanorama"}:
        meta["skip_integration_reason"] = "sample and condition are 1:1 collinear; integrating would erase the contrast"
    integrator = None if r_degraded else choose_integrator(meta, requested)
    if req_l in {"none", "off", "skip"}:
        meta["skip_integration_reason"] = "user disabled batch module"
    integrator_reason = explain_integrator(meta, requested, integrator)
    imputation = "none" if r_degraded else (state.get("imputation") or (cfg.get("modules") or {}).get("imputation") or "none")
    ambient = "none" if r_degraded else choose_ambient(
        meta.get("tissue"), state.get("ambient") or (cfg.get("modules") or {}).get("ambient")
    )
    intent = {} if r_degraded else parse_intent(state.get("user_query"), cfg)
    # Exploratory DEG intent alone must NOT hard-require pseudobulk (Reviewer would false-fail).
    wants_deg = False if r_degraded else bool(
        intent.get("condition_comparison") or "deg" in (intent.get("intents") or [])
    )
    if state.get("condition_key"):
        meta["condition_key"] = state["condition_key"]
        wants_deg = True
    if force_pseudobulk_de(meta):
        meta["force_pseudobulk_de"] = True
    # Hard gate only when condition_key + n_replicates≥2 (force_pseudobulk_de).
    needs_pb = False if r_degraded else bool(force_pseudobulk_de(meta))
    ct_model = None if r_degraded else choose_celltypist_model(meta.get("tissue"), meta.get("species"))
    resolution = state.get("resolution")
    if resolution is None:
        resolution = analysis_params(cfg).get("leiden_resolution")
    if r_degraded:
        route = ["plan_only"]
    else:
        intents = list(intent.get("intents") or ["qc", "clustering", "annotation"])
        if (wants_deg or needs_pb) and "deg" not in intents:
            intents.append("deg")
        from scagent.trajectory import should_plan_trajectory

        if should_plan_trajectory(meta, intent, cfg) and "trajectory" not in intents:
            intents.append("trajectory")
        if any(x in intents for x in ("clustering", "annotation", "deg")) and "enrichment" not in intents:
            intents.append("enrichment")
        route = resolve_route(
            intents,
            integrator=integrator,
            imputation=imputation,
            ambient=ambient,
            r_degraded=False,
        )

    skills = [] if r_degraded else recommend_skills(
        {
            **meta,
            "integrator": integrator,
            "user_query": state.get("user_query"),
            "task": state.get("user_query"),
            "intents": (intent or {}).get("intents") or [],
            "route": route,
        },
        language=language,
    )
    best_practices = [] if r_degraded else practices_for_route(
        route if not r_degraded else [],
        (intent or {}).get("intents") or [],
        state.get("user_query"),
    )
    rag = (
        ""
        if r_degraded
        else format_hits(
            retrieve_fused(
                f"{state.get('user_query') or ''} {meta.get('platform')} {meta.get('tissue')} "
                f"{' '.join(route)} scRNA-seq",
                route=route,
                intents=list((intent or {}).get("intents") or []),
                user_query=state.get("user_query"),
            )
        )
    )
    from agents.literature import fetch_phase_literature

    lit = (
        {"paper_excerpt": "", "paper_recs": []}
        if r_degraded
        else fetch_phase_literature(
            "plan",
            tissue=str(meta.get("tissue") or ""),
            platform=str(meta.get("platform") or ""),
            user_query=str(state.get("user_query") or ""),
        )
    )

    pref = parse_deg_preference(state.get("user_query"))
    tool_route = build_tool_route(
        meta,
        {"needs_pseudobulk": needs_pb, "integrator": integrator, "route": route if not r_degraded else []},
        cfg=cfg,
    )
    plan = {
        "objective": state.get("user_query") or "标准 scRNA-seq 分析",
        "species": meta.get("species"),
        "platform": meta.get("platform"),
        "n_samples": meta.get("n_samples"),
        "n_cells": meta.get("n_cells"),
        "need_batch_correction": bool(meta.get("need_batch_correction")),
        "integrator": integrator,
        "integrator_reason": integrator_reason,
        "sample_key": meta.get("sample_key"),
        "skip_integration_reason": meta.get("skip_integration_reason"),
        "imputation": imputation or "none",
        "ambient": ambient,
        "needs_pseudobulk": needs_pb,
        "force_pseudobulk_de": bool(meta.get("force_pseudobulk_de") or force_pseudobulk_de(meta)),
        "n_replicates": meta.get("n_replicates"),
        "deg_engine": resolve_forced_deg_engine(
            state.get("deg_engine") or pref.get("engine") or (cfg.get("deg") or {}).get("engine") or "auto"
        )
        if force_pseudobulk_de(meta)
        else (state.get("deg_engine") or pref.get("engine") or (cfg.get("deg") or {}).get("engine") or "auto"),
        "marker_method": state.get("marker_method") or pref.get("marker_method") or (cfg.get("deg") or {}).get("marker_method") or "auto",
        "deg_cross_validate": (
            state.get("deg_cross_validate")
            if state.get("deg_cross_validate") is not None
            else pref.get("cross_validate")
            if pref.get("cross_validate") is not None
            else (cfg.get("deg") or {}).get("cross_validate") or "auto"
        ),
        "condition_key": meta.get("condition_key") or state.get("condition_key"),
        "celltypist_model": ct_model,
        "intent": intent if not r_degraded else {"intents": [], "source": "r_degraded"},
        "resolution": resolution,
        "language": language,
        "r_degraded": r_degraded,
        "skills": skills,
        "best_practices": best_practices,
        "route": route,
        "dag": serialize_dag(route, integrator=integrator),
        "loop": "plan-and-solve",
        "collaboration": "multi-agent",
        "agents": assign_roles(route),
        "risks": [],
        "rag_excerpt": rag,
        "paper_excerpt": lit.get("paper_excerpt") or "",
        "paper_recs": lit.get("paper_recs") or [],
        "tool_route": tool_route,
    }
    if state.get("selection"):
        sel = state["selection"] or {}
        n_sel = sel.get("n") or len(sel.get("cell_ids") or [])
        plan["selection_n"] = n_sel
        plan["risks"].append(f"用户在交互 UMAP 框选了 {n_sel} 个细胞；针对该选区回答，不要当成全数据结论。")
    if r_degraded:
        plan["risks"].append(
            "现有 skills 为 Python/Scanpy。--language r 写出双重格式 Seurat Rmd，scAgent 不执行 R kernel。"
        )
    if meta.get("platform") == "parse":
        plan["risks"].append("Parse 平台 barcode 与 10x 不同，加载时不要套 Cell Ranger 默认假设。")
    if meta.get("batch_condition_confounded"):
        plan["risks"].append(
            "样本与条件 1:1 共线。auto 已跳过整合，避免把处理效应当批次抹掉（Luecken 2022 overcorrection）。"
        )
    if integrator:
        plan["risks"].append(
            f"多样本选用 {integrator}（n_cells={meta.get('n_cells')}, n_samples={meta.get('n_samples')}）。"
            "禁止把 UMAP 混匀当整合成功。"
        )
    if "trajectory" in (plan.get("route") or []):
        plan["risks"].append(
            "轨迹：先评估 PAGA 连续性；支持则拟合 DPT/Palantir 分化轴与基因趋势。"
            "scVelo 仅在 spliced/unspliced 时运行；Monocle3 需 R。推断轨迹不等于已验证的生物学命运。"
        )
    if plan.get("marker_method") == "mast":
        plan["risks"].append("MAST 是细胞水平 hurdle 模型，不是样本水平检验；组间结论仍走 pseudobulk + FDR。")
    if plan.get("deg_cross_validate") not in {False, "off", "false"}:
        plan["risks"].append("DEG/marker 将跑第二检验做交叉验证；共识基因比单方法列表更稳，仍不是因果。")
    if plan.get("force_pseudobulk_de"):
        plan["risks"].append(
            "检测到 condition_key 且生物学重复 n_replicates≥2：禁止 cell-level Wilcoxon 作组间结论，强制 pseudobulk + DESeq2/edgeR。"
        )
    if language != "r":
        llm = run_specialist(
            read_prompt("planner"),
            (
                f"用户任务: {state.get('user_query')}\n"
                f"metadata: {json.dumps(meta, ensure_ascii=False)}\n"
                f"integrator={integrator}\n"
                f"可选 skills（全部 bundled，* 为本次推荐）:\n{skill_catalog_text(metadata={**meta, 'integrator': integrator}, query=state.get('user_query'))}\n"
                f"knowledge/best_practices:\n{practices_catalog_text()}\n"
                f"本次选用 SOP: {', '.join(best_practices)}\n"
                f"RAG（SOP + 文献融合）:\n{rag}\n"
                f"文献段落（papers Methods/Results 加权）:\n{lit.get('paper_excerpt') or '（无）'}\n"
                "这是 Plan-and-Solve + 多智能体：Planner 只分工，QC / 聚类DEG / 解读 / 代码审计各司其职。"
                "禁止在 PCA/neighbors/UMAP/Leiden 之前做 DE 或 DPT/Monocle3。"
                "请结合文献最佳实践输出分析路线，并点名可引用的论文要点。"
            ),
        )
        plan["narrative"] = llm or (
            f"物种={plan['species']}，平台={plan['platform']}，样本数={plan['n_samples']}。"
            f"整合={integrator or '无'}。路线: {' → '.join(route)}。skills: {', '.join(skills)}。"
        )
    else:
        plan["narrative"] = (
            "R/Seurat：写出双重格式 analysis.Rmd（结论 + 可运行 Seurat 块）。"
            "scAgent 不执行 R kernel；仓库内 SOP 仍为 Scanpy。"
            " 建议路线: QC(MAD) → 归一化 → HVG → PCA → 整合(如需) → Leiden → 注释双验证。"
        )
    return plan
