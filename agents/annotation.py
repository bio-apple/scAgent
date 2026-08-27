from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from agents.markers import load_marker_catalog
from agents.templates import cluster_annotate_script
from rag.retriever import format_hits, retrieve


def build_annotation_plan(state: dict) -> dict:
    meta = dict(state.get("metadata") or {})
    if state.get("markers_path"):
        meta["markers_path"] = state["markers_path"]
    catalog = load_marker_catalog(meta.get("markers_path"), tissue=str(meta.get("tissue") or "pbmc"))
    rag = format_hits(
        retrieve(
            "CellTypist marker dual validation annotation",
            collections=["papers", "best_practices"],
        )
    )
    markers = format_hits(retrieve(str(meta.get("tissue") or "immune"), collection="markers"))
    code = cluster_annotate_script(meta, state.get("qc_strategy") or {}, state.get("plan") or {})
    plan = {
        "tiers": [
            "Leiden 无偏聚类",
            "参考映射（CellTypist 等）仅作假说",
            "层级 marker 双验证（≥2 阳性 + ≥1 阴性）",
            "冲突保留 marker；跨组织不让免疫模型覆盖",
        ],
        "forbid_single_gene": True,
        "low_confidence_cutoff": 0.5,
        "dual_validation": True,
        "catalog_tissue": catalog.get("tissue"),
        "n_cell_types": len(catalog.get("cell_types") or []),
        "rag_excerpt": rag,
        "marker_excerpt": markers,
        "code": code,
        "instructions": (
            "层级注释 Immune→T→CD8→Tex。"
            "自动注释（CellTypist/LLM）只是假说；≥2 阳性 + ≥1 阴性 marker 为生物学赋值。"
            "冲突保留 marker，标 annotation_conflict。跨组织禁用免疫模型覆盖。"
            "低置信且无 marker 时标 low_conf / unvalidated。"
        ),
        "hierarchical": True,
    }
    llm = run_specialist(
        read_prompt("annotation"),
        f"metadata={json.dumps(meta, ensure_ascii=False)}\nRAG:\n{rag}\nmarkers:\n{markers}",
    )
    if llm:
        plan["instructions"] = llm
    return plan
