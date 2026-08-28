"""Biological Interpretation Agent: GSEA/GSVA (or ORA) + RAG literature check."""

from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from rag.retriever import format_hits, retrieve_fused
from scagent.best_practices_loader import practices_for_phase
from scagent.kb import format_records, lookup_structured


def build_interpretation_plan(state: dict) -> dict:
    meta = state.get("metadata") or {}
    tissue = str(meta.get("tissue") or "default")
    plan = state.get("plan") or {}
    mets = (state.get("artifacts") or {}).get("metrics") or {}
    genes_hint = mets.get("deg_n_sig")
    rag = format_hits(
        retrieve_fused(
            f"{tissue} pathway enrichment GSEA GSVA Hallmark Heumos scRNA-seq",
            phase="interpret",
            route=list(plan.get("route") or []),
            user_query=state.get("user_query"),
            tissue=tissue,
        )
    )
    kb = format_records(
        lookup_structured(
            f"{tissue} pathway disease signature Hallmark",
            collections=["pathway", "disease_signature"],
            tissue=tissue,
            top_k=6,
        )
    )
    bp = practices_for_phase("interpret", route=list(plan.get("route") or []), query=state.get("user_query"))
    method = "gsea_or_ora"
    gsva_note = "GSVA 需 decoupler；未安装时只做 ORA/GSEA 风格的基因集过表达，不把缺失写成 GSVA 结果。"
    out = {
        "role": "bio_interpret",
        "agent": "Biological Interpretation Agent",
        "method": method,
        "gene_sets": "MSigDB Hallmark-like (offline subset); GSEA/GSVA if libraries present",
        "fdr_note": "GSEA 常用 FDR<0.25；ORA 报告 BH FDR。基因集选择比检验方法更关键（Heumos 2023）。",
        "gsva_note": gsva_note,
        "literature": rag,
        "knowledge": kb,
        "best_practices": bp,
        "deg_n_sig": genes_hint,
        "needs_de_input": True,
        "tissue": tissue,
        "integrator": plan.get("integrator"),
        "instructions": (
            "用 DEG 或 cluster marker 做通路富集；优先 GSEA，其次 ORA；GSVA 仅在有表达矩阵与 decoupler 时运行。"
            "文献验证走本地 RAG，不把 UMAP 混匀或单条通路 p 值写成机制结论。"
            "任何细胞状态断言（如耗竭 T）必须附证据链：≥2 marker（如 PDCD1+HAVCR2）+ GO/通路 p 值 + PubMed DOI/PMID。"
        ),
    }
    llm = run_specialist(
        read_prompt("bio_interpret"),
        (
            f"tissue={tissue}\n"
            f"deg_n_sig={genes_hint}\n"
            f"route={plan.get('route')}\n"
            f"RAG:\n{rag}\n"
            f"structured_kb:\n{kb}\n"
            "输出解读要点：用什么基因集、为何不是 GSVA 默认、如何对照文献。"
        ),
    )
    if llm:
        out["instructions"] = llm
    out["brief"] = json.dumps({k: out[k] for k in ("method", "gene_sets", "tissue") if k in out}, ensure_ascii=False)
    return out
