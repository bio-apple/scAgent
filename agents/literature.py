"""Phase-aware paper retrieval for Planner / QC / Annotation / Interpretation."""

from __future__ import annotations

from pathlib import Path

from rag.retriever import format_paper_hits, search_paper_knowledge

# Seed queries when user_query is short; still prepend tissue + user text.
PHASE_PAPER_QUERIES: dict[str, str] = {
    "plan": (
        "scRNA-seq analysis workflow best practices batch integration "
        "Harmony scVI quality control clustering"
    ),
    "qc": (
        "single-cell quality control mitochondrial proportion MAD threshold "
        "doublet detection ambient RNA filtering"
    ),
    "annotation": (
        "cell type annotation marker genes CellTypist reference mapping "
        "dual validation hierarchical labels"
    ),
    "interpret": (
        "pathway enrichment GSEA Hallmark differential expression "
        "pseudobulk biological interpretation cell state"
    ),
}


def _first_sentence(text: str, limit: int = 220) -> str:
    body = " ".join((text or "").split())
    if not body:
        return ""
    for sep in (". ", "; ", "。", "；"):
        if sep in body:
            body = body.split(sep, 1)[0] + ("." if sep.startswith(".") else "")
            break
    return body[:limit].rstrip(" ,;:") + ("…" if len(body) > limit else "")


def literature_recommendations(hits: list[dict], *, max_n: int = 4) -> list[str]:
    """Short bullets for reports: title — section — snippet."""
    out: list[str] = []
    for h in hits[:max_n]:
        stem = h.get("stem") or Path(str(h.get("source") or "")).stem
        sec = h.get("section") or "body"
        snip = _first_sentence(str(h.get("text") or ""))
        if not snip:
            continue
        out.append(f"{stem} [{sec}]: {snip}")
    return out


def fetch_phase_literature(
    phase: str,
    *,
    tissue: str = "",
    user_query: str = "",
    platform: str = "",
    top_k: int = 4,
) -> dict:
    """Return paper_excerpt (prompt), paper_recs (report bullets), and raw hits."""
    seed = PHASE_PAPER_QUERIES.get(phase) or PHASE_PAPER_QUERIES["plan"]
    query = " ".join(x for x in (user_query, tissue, platform, seed) if x).strip()
    try:
        hits = search_paper_knowledge(query, top_k=top_k)
    except Exception:
        hits = []
    return {
        "paper_hits": hits,
        "paper_excerpt": format_paper_hits(hits) if hits else "",
        "paper_recs": literature_recommendations(hits),
    }


def format_literature_report_block(
    *,
    plan: dict | None = None,
    qc: dict | None = None,
    ann: dict | None = None,
    interpret: dict | None = None,
    lang: str = "zh",
) -> str:
    """Markdown section for writer: literature-based best-practice suggestions."""
    zh = lang != "en"
    title = "## 文献最佳实践建议" if zh else "## Literature-based best practices"
    groups = [
        ("Planner", (plan or {}).get("paper_recs") or []),
        ("QC", (qc or {}).get("paper_recs") or []),
        ("Annotation", (ann or {}).get("paper_recs") or []),
        ("Interpretation", (interpret or {}).get("paper_recs") or []),
    ]
    lines = [title, ""]
    any_hit = False
    for name, recs in groups:
        if not recs:
            continue
        any_hit = True
        lines.append(f"### {name}")
        lines.append("")
        for r in recs:
            lines.append(f"- {r}")
        lines.append("")
    if not any_hit:
        lines.append(
            "- （未检索到文献段落。请先 `scagent parse-papers` 与 `scagent ingest`。）"
            if zh
            else "- (No paper passages retrieved. Run `scagent parse-papers` then `scagent ingest`.)"
        )
        lines.append("")
    else:
        lines.append(
            "以上摘录来自 `knowledge/papers/.parsed/`，用于支撑 QC / 注释 / 解读决策；"
            "不是对实验因果的证明。"
            if zh
            else "Excerpts from `knowledge/papers/.parsed/` support QC / annotation / interpretation decisions; "
            "they are not causal proof."
        )
        lines.append("")
    return "\n".join(lines)
