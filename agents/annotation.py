from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from agents.markers import choose_celltypist_model, load_marker_catalog
from agents.templates import cluster_annotate_script
from rag.retriever import format_hits, retrieve


def build_annotation_plan(state: dict) -> dict:
    meta = dict(state.get("metadata") or {})
    if state.get("markers_path"):
        meta["markers_path"] = state["markers_path"]
    catalog = load_marker_catalog(meta.get("markers_path"), tissue=str(meta.get("tissue") or "pbmc"))
    ct_model = choose_celltypist_model(meta.get("tissue"), meta.get("species"))
    plan_in = dict(state.get("plan") or {})
    plan_in.setdefault("celltypist_model", ct_model)
    rag = format_hits(
        retrieve(
            "CellTypist marker dual validation annotation",
            collections=["papers", "best_practices"],
        )
    )
    markers = format_hits(retrieve(str(meta.get("tissue") or "immune"), collection="markers"))
    code = cluster_annotate_script(meta, state.get("qc_strategy") or {}, plan_in)
    plan = {
        "tiers": [
            "Leiden 无偏聚类",
            f"参考映射（CellTypist {ct_model or '无匹配模型'}）仅作假说",
            "第二参考（SingleR/Azimuth/popV 或 marker Spearman）交叉验证",
            "层级 marker 双验证（≥2 阳性 + ≥1 阴性）",
            "冲突保留 marker；跨组织不让免疫模型覆盖",
        ],
        "forbid_single_gene": True,
        "low_confidence_cutoff": 0.5,
        "dual_validation": True,
        "celltypist_model": ct_model,
        "catalog_tissue": catalog.get("tissue"),
        "n_cell_types": len(catalog.get("cell_types") or []),
        "rag_excerpt": rag,
        "marker_excerpt": markers,
        "code": code,
        "instructions": (
            "层级注释 Immune→T→CD8→Tex。"
            "自动注释（CellTypist/SingleR）只是假说；≥2 阳性 + ≥1 阴性 marker 为生物学赋值。"
            "CellTypist 模型按组织选择；与第二参考交叉验证。冲突保留 marker。"
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
