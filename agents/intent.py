"""Structured analysis intent. LLM JSON when available; rule fallback never matches degradation as DEG."""

from __future__ import annotations

import json
import re
from typing import Any, Literal

from pydantic import BaseModel, Field

IntentName = Literal["qc", "clustering", "deg", "trajectory", "annotation", "enrichment"]
INTENT_NAMES: tuple[IntentName, ...] = ("qc", "clustering", "deg", "trajectory", "annotation", "enrichment")


class AnalysisIntent(BaseModel):
    intents: list[IntentName] = Field(default_factory=list)
    condition_comparison: bool = False
    source: str = "rules"


_DEG_POS = (
    r"差异表达",
    r"differential\s+express",
    r"\bpseudobulk\b",
    r"\bdeseq2?\b",
    r"\bedger\b",
    r"对照组",
    r"处理组",
    r"condition\s+vs",
    r"treated\s+vs",
    r"\bcase[- ]control\b",
    r"\bdeg\b",
    r"\bdge\b",
)
_DEG_NEG = (r"degrad", r"degree", r"\bedge\s+detect")
_TRAJ = (
    r"轨迹",
    r"伪时间",
    r"命运",
    r"分化轴",
    r"\bpaga\b",
    r"pseudotime",
    r"\btrajectory\b",
    r"diffusion\s+map",
    r"\bdpt\b",
    r"monocle",
    r"palantir",
    r"scvelo",
    r"\bvelocity\b",
    r"rna\s+velocity",
)
_ANN = (r"注释", r"cell\s*type", r"celltypist", r"annotat")
_CLUST = (r"聚类", r"\bleiden\b", r"\blouvain\b", r"cluster")
_QC = (r"质控", r"\bqc\b", r"doublet", r"mito", r"线粒体", r"scrublet")
_ENR = (r"通路", r"富集", r"\bgsea\b", r"\bgsva\b", r"enrichment", r"hallmark")


def _search(query: str, patterns: tuple[str, ...]) -> bool:
    return any(re.search(p, query, re.I) for p in patterns)


def rule_intent(query: str | None) -> dict[str, Any]:
    q = str(query or "")
    intents: list[str] = []
    if _search(q, _QC):
        intents.append("qc")
    if _search(q, _CLUST):
        intents.append("clustering")
    if _search(q, _ANN):
        intents.append("annotation")
    deg = _search(q, _DEG_POS) and not _search(q, _DEG_NEG)
    if deg:
        intents.append("deg")
    if _search(q, _TRAJ):
        intents.append("trajectory")
    if _search(q, _ENR):
        intents.append("enrichment")
    if not intents:
        intents = ["qc", "clustering", "annotation"]
    if any(x in intents for x in ("clustering", "annotation", "deg")) and "enrichment" not in intents:
        intents.append("enrichment")
    return {
        "intents": intents,
        "condition_comparison": bool(deg),
        "source": "rules",
    }


def parse_intent(query: str | None, cfg: dict | None = None) -> dict[str, Any]:
    """Prefer LLM JSON schema; fall back to rules (no 'deg' substring trap)."""
    rules = rule_intent(query)
    llm = _llm_intent(query, cfg)
    if not llm:
        return rules
    names = [x for x in (llm.get("intents") or []) if x in INTENT_NAMES]
    if not names:
        names = rules["intents"]
    cond = bool(llm.get("condition_comparison"))
    if "deg" in names:
        cond = True
    if _search(str(query or ""), _DEG_NEG):
        names = [x for x in names if x != "deg"]
        cond = False
    return {"intents": names, "condition_comparison": cond, "source": "llm"}


def _llm_intent(query: str | None, cfg: dict | None) -> dict[str, Any] | None:
    try:
        from agents.common import invoke_json, read_prompt
    except Exception:
        return None
    from agents.common import get_llm

    if get_llm(cfg) is None:
        return None
    system = (
        read_prompt("planner")
        + "\n\nReturn ONLY a JSON object with keys:\n"
        '{"intents": ["qc"|"clustering"|"deg"|"trajectory"|"annotation"|"enrichment", ...], '
        '"condition_comparison": boolean}\n'
        "Do not set deg for RNA degradation, degree, or edge detection. "
        "Standard annotation/QC tasks should not include deg unless a condition contrast is requested."
    )
    data = invoke_json(system, f"用户任务: {query}", cfg=cfg)
    if not isinstance(data, dict):
        return None
    try:
        return AnalysisIntent.model_validate(data).model_dump()
    except Exception:
        return None
