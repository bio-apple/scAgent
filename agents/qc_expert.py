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
    meta = state.get("metadata") or {}
    tissue = str(meta.get("tissue") or "default")
    prof = _profile(tissue)
    rag = format_hits(retrieve(f"{tissue} mitochondrial QC MAD violin scatter", collection="papers"))
    strategy = {
        "tissue": tissue,
        "nmads": 5 if tissue not in {"heart", "tumor"} else 6,
        "plots_required": ["violin", "scatter", "mad"],
        "metrics": ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        "extra_qc": list(prof.get("extra_qc") or []),
        "pct_mt_note": prof.get("pct_mt_note", ""),
        "doublets": True,
        "rag_excerpt": rag,
        "protocol": (
            "1. calculate_qc_metrics（mt/ribo/hb）\n"
            "2. Violin：n_genes_by_counts, total_counts, pct_counts_mt\n"
            "3. Scatter：counts vs genes；counts vs pct_mt\n"
            "4. MAD 离群（默认 5 MAD；心脏/肿瘤放宽）并记录移除数\n"
            "5. 双细胞检测后再注释"
        ),
    }
    llm = run_specialist(
        read_prompt("qc_expert"),
        f"metadata={json.dumps(meta, ensure_ascii=False)}\nprofile={json.dumps(prof, ensure_ascii=False)}\nRAG:\n{rag}",
    )
    if llm:
        strategy["protocol"] = llm
    return strategy
