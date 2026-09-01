from __future__ import annotations

import json

from agents.common import read_prompt, run_specialist
from rag.retriever import format_hits, retrieve_fused
from scagent.best_practices_loader import practices_for_phase
from scagent.config import load_config
from scagent.doublets import resolve_doublet_methods
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
    override = dict(state.get("qc_override") or {})
    if state.get("qc_choice") and state.get("hitl_mt"):
        from scagent.hitl import mt_override

        override.update(mt_override(state["hitl_mt"], state.get("qc_choice")))
    elif state.get("qc_choice"):
        from scagent.hitl import build_mt_decision, mt_override

        override.update(mt_override(build_mt_decision(state), state.get("qc_choice")))
    prof.update(override)
    if state.get("qc_method"):
        qc_cfg["method"] = state["qc_method"]
    rag = format_hits(
        retrieve_fused(
            f"{tissue} mitochondrial QC MAD percentile violin scatter",
            phase="qc",
            user_query=state.get("user_query"),
        )
    )
    from agents.literature import fetch_phase_literature

    lit = fetch_phase_literature(
        "qc",
        tissue=tissue,
        user_query=str(state.get("user_query") or ""),
    )
    bp = practices_for_phase("qc", query=state.get("user_query"))
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
    doublet_filter = str(state.get("doublet_filter") or qc_cfg_full.get("doublet_filter") or "high_conf")
    regress_cc = state.get("regress_cell_cycle") or qc_cfg_full.get("regress_cell_cycle") or "auto"
    normalization = str(
        state.get("normalization") or qc_cfg_full.get("normalization") or "log1p"
    ).lower()
    per_sample_qc = qc_cfg_full.get("per_sample")
    if state.get("per_sample_qc") is not None:
        per_sample_qc = bool(state.get("per_sample_qc"))
    extra_qc = list(prof.get("extra_qc") or [])
    if meta.get("need_hb_qc") and "hb" not in extra_qc:
        extra_qc.append("hb")
    warn_pct = int(qc_cfg.get("overfilter_warn_pct") or 30)
    n_samples = int(meta.get("n_samples") or 1)
    if meta.get("need_batch_correction"):
        n_samples = max(n_samples, 2)
    doublet_req = str(state.get("doublet_methods") or qc_cfg.get("doublet_methods") or "auto")
    doublet_resolved = resolve_doublet_methods(doublet_req, tissue=tissue, n_samples=n_samples)
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
        "doublet_filter": doublet_filter,
        "doublet_methods": doublet_req,
        "doublet_methods_resolved": doublet_resolved,
        "ambient": ambient,
        "ambient_requested": ambient_req,
        "normalization": normalization,
        "per_sample_qc": per_sample_qc,
        "regress_cell_cycle": regress_cc,
        "mt_mad_side": "high",
        "overfilter_warn_pct": warn_pct,
        "imputation": imputation,
        "hitl_choice": state.get("qc_choice"),
        "rag_excerpt": rag,
        "paper_excerpt": lit.get("paper_excerpt") or "",
        "paper_recs": lit.get("paper_recs") or [],
        "best_practices": bp,
        "protocol": (
            f"QC method={method}（组织={tissue}）。禁止默认 mito%<5 / nFeature>200。\n"
            "1. calculate_qc_metrics(log1p=True) + Violin/Scatter\n"
            "2. MAD 和/或 percentile 动态阈值（多样本时 per-sample）；硬阈值仅当 config.qc.hard 非 null\n"
            f"3. pct_mt 单侧高；记录过滤比例，>{warn_pct}% 警告；组织说明：{prof.get('pct_mt_note') or ''}\n"
            f"4. 双细胞：Scrublet"
            + (
                " + scDblFinder（无 R 则 count-simulation 启发式，非真 scDblFinder）→ doublet_call 三级："
                "high_conf（两法一致且 score>0.8）| low_conf（仅一法）| singlet"
                if len(doublet_resolved) > 1
                else " → doublet_call：high_conf（score>0.8）| low_conf | singlet"
            )
            + f"；remove_doublets={remove_doublets} doublet_filter={doublet_filter}\n"
            f"5. ambient={ambient}（auto 在无真 SoupX 时为 none；显式 soupx 不可用则不改 counts）"
            f"；normalization={normalization}；cell cycle regress={regress_cc}\n"
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
        (
            f"metadata={json.dumps(meta, ensure_ascii=False)}\n"
            f"profile={json.dumps(prof, ensure_ascii=False)}\n"
            f"RAG:\n{rag}\n"
            f"文献段落:\n{lit.get('paper_excerpt') or '（无）'}\n"
            "请输出 QC 协议，并引用文献中的 mito%/MAD/双细胞相关建议。"
        ),
    )
    if llm:
        strategy["protocol"] = llm
    return strategy
