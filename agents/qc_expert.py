from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from rag.retriever import format_hits, retrieve
from scagent.config import load_config


def _profile(tissue: str) -> dict:
    cfg = load_config()
    profiles = cfg.get("qc_profiles") or {}
    key = (tissue or "default").lower()
    return profiles.get(key) or profiles.get("default") or {}


def build_qc_strategy(state: dict) -> dict:
    cfg = load_config()
    qc_cfg = dict(cfg.get("qc") or {})
    meta = state.get("metadata") or {}
    tissue = str(meta.get("tissue") or "default")
    prof = dict(_profile(tissue))
    override = state.get("qc_override") or {}
    prof.update(override)
    if state.get("qc_method"):
        qc_cfg["method"] = state["qc_method"]
    rag = format_hits(
        retrieve(
            f"{tissue} mitochondrial QC MAD percentile violin scatter",
            collections=["papers", "best_practices"],
        )
    )
    nmads = int(prof.get("nmads") or (6 if tissue in {"heart", "tumor", "kidney"} else 5))
    method = str(qc_cfg.get("method") or "mad")
    modules = cfg.get("modules") or {}
    imputation = state.get("imputation") or modules.get("imputation") or "none"
    strategy = {
        "tissue": tissue,
        "method": method,
        "nmads": nmads,
        "percentile": dict(qc_cfg.get("percentile") or {}),
        "hard": dict(qc_cfg.get("hard") or {}),
        "plots_required": ["violin", "scatter", "mad"],
        "metrics": ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        "extra_qc": list(prof.get("extra_qc") or []),
        "pct_mt_note": prof.get("pct_mt_note", ""),
        "doublets": True,
        "ambient": tissue in {"brain", "tumor"} or bool(prof.get("ambient")),
        "mt_mad_side": "high",
        "overfilter_warn_pct": int(qc_cfg.get("overfilter_warn_pct") or 30),
        "imputation": imputation,
        "rag_excerpt": rag,
        "protocol": (
            f"QC method={method}（组织={tissue}）。禁止默认 mito%<5 / nFeature>200。\n"
            "1. calculate_qc_metrics(log1p=True) + Violin/Scatter\n"
            "2. MAD 和/或 percentile 动态阈值；硬阈值仅当 config.qc.hard 非 null\n"
            "3. pct_mt 单侧高；记录过滤比例，>30% 警告\n"
            f"4. Scrublet；imputation={imputation}（写入 layers['imputed']，不覆盖 X）"
        ),
    }
    llm = run_specialist(
        read_prompt("qc_expert"),
        f"metadata={json.dumps(meta, ensure_ascii=False)}\nprofile={json.dumps(prof, ensure_ascii=False)}\nRAG:\n{rag}",
    )
    if llm:
        strategy["protocol"] = llm
    return strategy
