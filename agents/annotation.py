from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from rag.retriever import format_hits, retrieve


def build_annotation_plan(state: dict) -> dict:
    meta = state.get("metadata") or {}
    rag = format_hits(retrieve("CellTypist marker dual validation annotation", collection="papers"))
    markers = format_hits(retrieve(str(meta.get("tissue") or "immune"), collection="markers"))
    plan = {
        "tiers": [
            "Leiden 无偏聚类",
            "CellTypist 参考映射 + majority vote + 置信度",
            "至少两个独立 marker + 一个阴性 marker",
        ],
        "forbid_single_gene": True,
        "low_confidence_cutoff": 0.5,
        "rag_excerpt": rag,
        "marker_excerpt": markers,
        "instructions": (
            "按 single-cell-annotation-guide 三层策略注释。"
            "CellTypist 之后必须用 canonical marker 双向验证。"
            "禁止单基因定论。"
        ),
    }
    llm = run_specialist(
        read_prompt("annotation"),
        f"metadata={json.dumps(meta, ensure_ascii=False)}\nRAG:\n{rag}\nmarkers:\n{markers}",
    )
    if llm:
        plan["instructions"] = llm
    return plan
