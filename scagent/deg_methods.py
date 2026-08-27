"""DEG method preference and gene-list overlap (cluster markers + pseudobulk)."""

from __future__ import annotations

import re
from typing import Any, Iterable

_SCANPY_METHODS = {"wilcoxon": "wilcoxon", "t-test": "t-test", "ttest": "t-test", "t_test": "t-test", "t test": "t-test"}
_ENGINES = ("edger", "deseq2", "ttest")


def parse_deg_preference(query: str | None) -> dict[str, Any]:
    """Read method preference from the user query. Does not by itself imply condition DE."""
    q = str(query or "")
    low = q.lower()
    engine = None
    marker = None
    cv: bool | None = None
    if re.search(r"\bdeseq2?\b", low):
        engine = "deseq2"
    elif re.search(r"\bedger\b", low):
        engine = "edger"
    if re.search(r"\bmast\b", low):
        marker = "mast"
    if re.search(r"t\s*检验|t[-_ ]?test|\bttest\b", low):
        marker = marker or "t-test"
        if re.search(r"pseudobulk|组间|condition", low):
            engine = engine or "ttest"
    if re.search(r"wilcox|秩和|mann[- ]?whitney", low):
        marker = marker or "wilcoxon"
    if re.search(r"只用|仅用|only\s+(use\s+)?", q, re.I) and (engine or marker):
        cv = False
    if re.search(r"交叉验证|cross[-\s]?valid|多种(检验|方法)|共识基因", q, re.I):
        cv = True
    if marker == "mast" and re.search(r"wilcox|t[-_ ]?test|deseq|edger", low):
        cv = True if cv is None else cv
    if re.search(r"wilcox", low) and re.search(r"t[-_ ]?test|t检验", low):
        cv = True if cv is None else cv
    return {"engine": engine, "marker_method": marker, "cross_validate": cv}


def resolve_marker_method(raw: str | None) -> str:
    s = str(raw or "auto").lower().strip()
    if s in {"auto", "", "none", "default"}:
        return "wilcoxon"
    if s == "mast":
        return "mast"
    return _SCANPY_METHODS.get(s, "wilcoxon")


def resolve_cross_validate(raw: Any, *, explicit: bool | None = None) -> bool:
    if explicit is True:
        return True
    if explicit is False:
        return False
    s = str(raw if raw is not None else "auto").lower().strip()
    if s in {"off", "false", "0", "no", "none"}:
        return False
    if s in {"on", "true", "1", "yes", "always", "both"}:
        return True
    return True


def alt_scanpy_method(primary: str) -> str:
    return "t-test" if primary == "wilcoxon" else "wilcoxon"


def alt_engine(primary: str) -> str:
    p = str(primary or "auto").lower()
    if "deseq" in p:
        return "edger"
    if "ttest" in p or p == "t-test":
        return "edger"
    return "deseq2"


def gene_overlap(a: Iterable[str], b: Iterable[str]) -> dict[str, Any]:
    sa = {str(x) for x in a if x}
    sb = {str(x) for x in b if x}
    inter = sa & sb
    union = sa | sb
    return {
        "n_a": len(sa),
        "n_b": len(sb),
        "n_overlap": len(inter),
        "jaccard": round((len(inter) / len(union)) if union else 0.0, 4),
        "overlap_genes": sorted(inter)[:200],
    }


def sig_genes(rows: list[dict], *, fdr: float = 0.05, gene_key: str = "gene") -> list[str]:
    out = []
    for r in rows or []:
        try:
            q = float(r.get("fdr", 1))
        except (TypeError, ValueError):
            continue
        if q <= fdr and r.get(gene_key):
            out.append(str(r[gene_key]))
    return out
