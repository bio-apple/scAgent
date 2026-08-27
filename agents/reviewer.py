from __future__ import annotations

import re
from pathlib import Path

from agents.common import read_prompt, run_specialist
from agents.templates import LOCKED_END, LOCKED_START


def audit_code(code: str, metadata: dict | None = None, phase: str = "qc") -> dict:
    """Deterministic statistical audit. LLM cannot override hard fails."""
    metadata = metadata or {}
    text = code or ""
    low = text.lower()
    issues: list[str] = []

    has_violin = "violin" in low or "vlnplot" in low
    has_scatter = "scatter" in low
    has_mad = "median_abs_deviation" in low or bool(re.search(r"\bmad\b", low))
    has_locked = LOCKED_START in text and LOCKED_END in text
    mt_one_sided = 'side="high"' in text or "side='high'" in text or 'side = "high"' in text
    log1p_qc = "log1p=true" in low or "log1p = true" in low
    has_celltypist = "celltypist" in low
    has_dual = ("positive" in low and "negative" in low) or "dual" in low
    has_seed = "seed" in low
    has_pseudobulk_note = "pseudobulk" in low and ("fdr" in low or "padj" in low or "多重" in text)

    if phase == "qc":
        if not has_violin:
            issues.append("QC 缺少 Violin/VlnPlot")
        if not has_scatter:
            issues.append("QC 缺少 Scatter")
        if not has_mad:
            issues.append("QC 未使用 MAD 自适应阈值")
        if not has_locked:
            issues.append("QC 缺少不可删除的 LOCKED QC 代码块")
        if has_mad and not mt_one_sided:
            issues.append("线粒体 MAD 应为单侧（只滤高 pct_mt）")
        if not log1p_qc:
            issues.append("calculate_qc_metrics 需 log1p=True，以便 MAD 使用 log1p_total_counts 列")
        if re.search(r"pct[_]?mt\s*[<>=]+\s*(5|10)\b", low):
            issues.append("使用了固定 pctMT=5/10 阈值，应改为看分布 + MAD/percentile")

    if phase == "downstream":
        if metadata.get("need_batch_correction"):
            integrated = any(k in low for k in ("harmony", "scvi", "scanorama", "cca"))
            skipped = "skip integration" in low or "跳过整合" in text or "user disabled batch" in low
            if not integrated and not skipped:
                issues.append("多样本未做整合，也未声明跳过理由")
        if "umap" in low and re.search(r"leiden\(.*umap|cluster.*umap", low):
            issues.append("疑似在 UMAP 坐标上聚类")
        if not has_celltypist:
            issues.append("注释未调用 CellTypist（或未显式失败降级）")
        if not has_dual:
            issues.append("注释缺少 dual validation（≥2 阳性 + ≥1 阴性 marker）")
        if "cell_type_l1" not in low and "lineage" not in low:
            issues.append("注释缺少层级字段（cell_type_l1 / lineage）")
        if "rank_genes_groups" in low and not has_pseudobulk_note:
            issues.append("差异表达未声明探索-only / 组间须 pseudobulk+FDR")
        if "cell_type" in low and not has_dual:
            issues.append("单基因或无双验证的细胞类型赋值")

    if not has_seed:
        issues.append("未固定随机种子")

    hard = [x for x in issues if not x.startswith("未固定")]
    if phase == "qc":
        passed = bool(has_violin and has_scatter and has_mad and has_locked) and not any(
            x.startswith("QC ") or "pctMT=10" in x or "单侧" in x or "log1p=True" in x for x in issues
        )
    else:
        passed = not any(
            x.startswith("注释") or "多样本未做整合" in x or "UMAP 坐标" in x or "pseudobulk" in x or "层级" in x
            for x in issues
        )
        if not (has_celltypist and has_dual):
            passed = False

    return {
        "passed": passed,
        "issues": issues,
        "required_fixes": hard if not passed else [],
        "has_violin": has_violin,
        "has_scatter": has_scatter,
        "has_mad": has_mad,
        "has_locked": has_locked,
        "has_celltypist": has_celltypist,
        "has_dual": has_dual,
        "phase": phase,
    }


def audit_execution(
    execution: dict | None,
    artifacts: dict | None,
    *,
    phase: str,
    execute_code: bool,
    metadata: dict | None = None,
) -> dict:
    execution = execution or {}
    artifacts = artifacts or {}
    metadata = metadata or {}
    issues: list[str] = []
    if not execute_code or not execution.get("executed"):
        return {
            "passed": True,
            "skipped": True,
            "issues": ["未执行代码，结果指标未验证"],
        }

    if not execution.get("ok"):
        issues.append("执行失败（非零 returncode）")
        err = (execution.get("stderr") or "")[-800:]
        if err:
            issues.append(f"stderr: {err[:300]}")

    metrics = artifacts.get("metrics") or {}
    # phase-specific metrics may be nested
    phase_art = (artifacts.get("phases") or {}).get(phase) or {}
    metrics = {**metrics, **(phase_art.get("metrics") or {})}

    if phase == "qc":
        pct = metrics.get("pct_removed")
        if pct is not None and float(pct) > 30:
            issues.append(f"过度过滤：MAD 移除 {pct:.1f}% 细胞")
        h5ad = (artifacts.get("h5ads") or {}).get("qc") or (phase_art.get("h5ads") or {}).get("qc")
        if execution.get("ok") and not h5ad:
            issues.append("未写出 adata_qc.h5ad")
        figs = artifacts.get("figures") or phase_art.get("figures") or []
        names = " ".join(Path(p).name.lower() for p in figs)
        if execution.get("ok"):
            if "violin" not in names:
                issues.append("执行后缺少 violin 图")
            if "scatter" not in names:
                issues.append("执行后缺少 scatter 图")

    if phase == "downstream":
        h5ad = (artifacts.get("h5ads") or {}).get("processed")
        if execution.get("ok") and not h5ad:
            issues.append("未写出 adata_processed.h5ad")
        mix = metrics.get("batch_cluster_dominance")
        if metadata.get("need_batch_correction") and mix is not None and float(mix) >= 0.95:
            issues.append(f"整合质量可疑：cluster 内主导批次比例 {mix:.2f}（可能未混合）")

    passed = execution.get("ok") is True and not any(
        x.startswith("执行失败") or "缺少" in x or "未写出" in x for x in issues
    )
    # overfilter is warning: still fail so coder/user sees it, but we treat as fail
    if any("过度过滤" in x for x in issues):
        passed = False
    return {"passed": passed, "skipped": False, "issues": issues, "metrics": metrics}


def review_state(state: dict, phase: str | None = None) -> dict:
    phase = phase or state.get("phase") or "qc"
    code = state.get("code") or (state.get("code_qc") if phase == "qc" else state.get("code_downstream")) or ""
    meta = state.get("metadata") or {}
    code_result = audit_code(code, meta, phase=phase)
    exe = state.get("execution") or (
        state.get("execution_qc") if phase == "qc" else state.get("execution_downstream")
    )
    exe_result = audit_execution(
        exe,
        state.get("artifacts"),
        phase=phase,
        execute_code=bool(state.get("execute_code")),
        metadata=meta,
    )
    issues = list(code_result.get("issues") or []) + list(exe_result.get("issues") or [])
    passed = bool(code_result.get("passed")) and bool(exe_result.get("passed"))
    result = {
        **code_result,
        "passed": passed,
        "issues": issues,
        "required_fixes": issues if not passed else [],
        "execution_audit": exe_result,
        "phase": phase,
    }
    llm = run_specialist(
        read_prompt("reviewer"),
        f"phase={phase}\ncode:\n{code[:8000]}\nmetadata={meta}\naudit={result}",
    )
    if llm:
        result["narrative"] = llm
    return result
