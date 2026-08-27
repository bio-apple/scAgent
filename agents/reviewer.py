from __future__ import annotations

import re

from agents.common import read_prompt, run_specialist


def audit_code(code: str, metadata: dict | None = None) -> dict:
    """Deterministic statistical audit. LLM may add comments but cannot override hard fails."""
    metadata = metadata or {}
    text = code or ""
    low = text.lower()
    issues: list[str] = []

    has_violin = "violin" in low or "vlnplot" in low
    has_scatter = "scatter" in low
    has_mad = "median_abs_deviation" in low or bool(re.search(r"\bmad\b", low))
    if not has_violin:
        issues.append("QC 缺少 Violin/VlnPlot")
    if not has_scatter:
        issues.append("QC 缺少 Scatter")
    if not has_mad:
        issues.append("QC 未使用 MAD 自适应阈值")

    if re.search(r"pct[_]?mt\s*[<>=]+\s*10", low):
        issues.append("使用了固定 pctMT=10 阈值，应改为看分布 + MAD")

    if metadata.get("need_batch_correction"):
        if "harmony" not in low and "scvi" not in low and "integrate" not in low:
            issues.append("多样本未做整合，也未声明跳过理由")

    if "umap" in low and re.search(r"leiden\(.*umap|cluster.*umap", low):
        issues.append("疑似在 UMAP 坐标上聚类")

    if "rank_genes_groups" in low and "pseudobulk" not in low and "fdr" not in low and "padj" not in low:
        # Wilcoxon 探索可以，但必须声明；缺声明则警告而非直接 fail
        issues.append("差异表达未声明多重校正/pseudobulk 用途（探索 vs 组间结论）")

    passed = not any(
        x.startswith("QC ") or "单基因" in x or "多样本未做整合" in x for x in issues
    )
    # QC trio is mandatory
    if not (has_violin and has_scatter and has_mad):
        passed = False

    return {
        "passed": passed,
        "issues": issues,
        "required_fixes": issues if not passed else [],
        "has_violin": has_violin,
        "has_scatter": has_scatter,
        "has_mad": has_mad,
    }


def review_state(state: dict) -> dict:
    code = state.get("code") or ""
    meta = state.get("metadata") or {}
    result = audit_code(code, meta)
    llm = run_specialist(
        read_prompt("reviewer"),
        f"code:\n{code[:8000]}\nmetadata={meta}\n deterministic_audit={result}",
    )
    if llm:
        result["narrative"] = llm
    return result
