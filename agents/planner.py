from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from rag.retriever import format_hits, retrieve
from scagent.skills_loader import recommend_skills, skill_catalog_text


def build_plan(state: dict) -> dict:
    meta = state.get("metadata") or {}
    language = state.get("language") or "python"
    skills = recommend_skills(meta, language=language)
    rag = format_hits(
        retrieve(
            f"{meta.get('platform')} {meta.get('tissue')} scRNA-seq best practices integration",
            collection="papers",
        )
    )
    route = ["qc", "normalize", "hvg", "pca", "neighbors", "leiden", "umap", "annotate"]
    if meta.get("need_batch_correction"):
        route.insert(route.index("neighbors"), "harmony")
    plan = {
        "objective": state.get("user_query") or "标准 scRNA-seq 分析",
        "species": meta.get("species"),
        "platform": meta.get("platform"),
        "n_samples": meta.get("n_samples"),
        "need_batch_correction": bool(meta.get("need_batch_correction")),
        "language": language,
        "skills": skills,
        "route": route,
        "risks": [],
        "rag_excerpt": rag,
    }
    if meta.get("platform") == "parse":
        plan["risks"].append("Parse 平台 barcode 与 10x 不同，加载时不要套 Cell Ranger 默认假设。")
    if meta.get("need_batch_correction"):
        plan["risks"].append("多样本：整合前确认批次是技术噪声还是生物学；禁止把 UMAP 混匀当成功。")
    if language == "r":
        plan["risks"].append("仓库现有 skills 以 Python/Scanpy 为主；R 路径缺少同等 SOP，需人工核验。")

    llm = run_specialist(
        read_prompt("planner"),
        (
            f"用户任务: {state.get('user_query')}\n"
            f"metadata: {json.dumps(meta, ensure_ascii=False)}\n"
            f"可选 skills:\n{skill_catalog_text()}\n"
            f"RAG:\n{rag}\n"
            "请输出分析路线。"
        ),
    )
    if llm:
        plan["narrative"] = llm
    else:
        plan["narrative"] = (
            f"物种={plan['species']}，平台={plan['platform']}，样本数={plan['n_samples']}。"
            f"路线: {' → '.join(route)}。启用 skills: {', '.join(skills)}。"
        )
    return plan
