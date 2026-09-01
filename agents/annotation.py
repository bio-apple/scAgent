from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from agents.markers import choose_celltypist_model, load_marker_catalog
from agents.templates import cluster_annotate_script
from rag.retriever import format_hits, retrieve_fused
from scagent.best_practices_loader import practices_for_phase
from scagent.kb import format_records, lookup_structured


def build_annotation_plan(state: dict) -> dict:
    meta = dict(state.get("metadata") or {})
    if state.get("markers_path"):
        meta["markers_path"] = state["markers_path"]
    catalog = load_marker_catalog(meta.get("markers_path"), tissue=str(meta.get("tissue") or "unknown"))
    ct_model = choose_celltypist_model(meta.get("tissue"), meta.get("species"))
    plan_in = dict(state.get("plan") or {})
    if state.get("resolution") is not None:
        plan_in["resolution"] = state["resolution"]
    plan_in.setdefault("celltypist_model", ct_model)
    tissue = str(meta.get("tissue") or "unknown")
    rag = format_hits(
        retrieve_fused(
            f"{tissue} CellTypist marker dual validation annotation",
            phase="annotation",
            route=list(plan_in.get("route") or []),
            user_query=state.get("user_query"),
            tissue=tissue,
            include_markers=True,
        )
    )
    from agents.literature import fetch_phase_literature

    lit = fetch_phase_literature(
        "annotation",
        tissue=tissue,
        user_query=str(state.get("user_query") or ""),
    )
    bp = practices_for_phase("annotation", route=list(plan_in.get("route") or []), query=state.get("user_query"))
    markers = format_records(
        lookup_structured(
            f"{tissue} cell type markers ontology",
            collections=["marker_db", "cell_ontology", "tissue_reference"],
            tissue=tissue,
            top_k=8,
        )
    )
    code = cluster_annotate_script(meta, state.get("qc_strategy") or {}, plan_in)
    plan = {
        "tiers": [
            "Leiden 无偏聚类",
            f"参考映射（CellTypist {ct_model or '无匹配模型'}）+ scANVI 后备（max_prob<0.8 → scagent_annotation）",
            "独立证据：cluster Wilcoxon 基因 ∩ catalog；可选 SingleR/popV",
            "层级 marker 双验证（≥2 阳性 + ≥1 阴性）",
            "三路表决融合：≥2 一致才赋值；冲突标 mixed；单路 auto 标 unvalidated",
        ],
        "forbid_single_gene": True,
        "low_confidence_cutoff": 0.8,
        "dual_validation": True,
        "celltypist_model": ct_model,
        "catalog_tissue": catalog.get("tissue"),
        "catalog_warning": catalog.get("warning"),
        "n_cell_types": len(catalog.get("cell_types") or []),
        "rag_excerpt": rag,
        "paper_excerpt": lit.get("paper_excerpt") or "",
        "paper_recs": lit.get("paper_recs") or [],
        "best_practices": bp,
        "marker_excerpt": markers,
        "code": code,
        "instructions": (
            "自动注释：CellTypist 先行；max_prob<0.8 的细胞触发 scANVI 半监督，写入 obs['scagent_annotation']。"
            "独立证据：cluster DE∩catalog + ≥2 阳性/≥1 阴性 marker。"
            "fuse_annotation 多数表决：marker + scagent_annotation + deg_label；≥2 路一致才定 cell_type。"
        ),
        "hierarchical": True,
    }
    llm = run_specialist(
        read_prompt("annotation"),
        (
            f"metadata={json.dumps(meta, ensure_ascii=False)}\n"
            f"RAG:\n{rag}\n"
            f"文献段落:\n{lit.get('paper_excerpt') or '（无）'}\n"
            f"markers:\n{markers}\n"
            "请输出注释指令，并引用文献中的 marker / 参考映射最佳实践。"
        ),
    )
    if llm:
        plan["instructions"] = llm
    return plan
