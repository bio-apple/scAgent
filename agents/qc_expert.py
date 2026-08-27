from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from rag.retriever import format_hits, retrieve
from scagent.config import load_config
from scagent.preprocess import choose_ambient


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
    ambient_req = state.get("ambient") or modules.get("ambient") or "auto"
    ambient = choose_ambient(tissue, ambient_req)
    qc_cfg_full = cfg.get("qc") or {}
    remove_doublets = bool(state.get("remove_doublets"))
    if state.get("remove_doublets") is None:
        remove_doublets = bool(qc_cfg_full.get("remove_doublets"))
    regress_cc = state.get("regress_cell_cycle") or qc_cfg_full.get("regress_cell_cycle") or "auto"
    extra_qc = list(prof.get("extra_qc") or [])
    if meta.get("need_hb_qc") and "hb" not in extra_qc:
        extra_qc.append("hb")
    warn_pct = int(qc_cfg.get("overfilter_warn_pct") or 30)
    strategy = {
        "tissue": tissue,
        "method": method,
        "nmads": nmads,
        "percentile": dict(qc_cfg.get("percentile") or {}),
        "hard": dict(qc_cfg.get("hard") or {}),
        "plots_required": ["violin", "scatter", "mad"],
        "metrics": ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        "extra_qc": extra_qc,
        "pct_mt_note": prof.get("pct_mt_note", ""),
        "doublets": True,
        "remove_doublets": remove_doublets,
        "ambient": ambient,
        "ambient_requested": ambient_req,
        "regress_cell_cycle": regress_cc,
        "mt_mad_side": "high",
        "overfilter_warn_pct": warn_pct,
        "imputation": imputation,
        "rag_excerpt": rag,
        "protocol": (
            f"QC method={method}（组织={tissue}）。禁止默认 mito%<5 / nFeature>200。\n"
            "1. calculate_qc_metrics(log1p=True) + Violin/Scatter\n"
            "2. MAD 和/或 percentile 动态阈值；硬阈值仅当 config.qc.hard 非 null\n"
            f"3. pct_mt 单侧高；记录过滤比例，>{warn_pct}% 警告；组织说明：{prof.get('pct_mt_note') or ''}\n"
            f"4. Scrublet → predicted_doublet；remove_doublets={remove_doublets}\n"
            f"5. ambient={ambient}；cell cycle score + regress={regress_cc}\n"
            f"6. imputation={imputation}（写入 layers['imputed']，不覆盖 X）"
        ),
    }
    layer = meta.get("expression_layer")
    if layer and layer != "counts":
        strategy["expression_layer"] = layer
        strategy["protocol"] += (
            f"\n7. 输入 X 检测为 {layer}：禁止重复 normalize/scale。"
            " log1p 则跳过归一化；scaled 需从 layers['counts'] 恢复。"
        )
    llm = run_specialist(
        read_prompt("qc_expert"),
        f"metadata={json.dumps(meta, ensure_ascii=False)}\nprofile={json.dumps(prof, ensure_ascii=False)}\nRAG:\n{rag}",
    )
    if llm:
        strategy["protocol"] = llm
    return strategy
