from __future__ import annotations

from datetime import date

from agents.reviewer import format_review_card, publication_review


CAPTION_EN = {
    "violin": "QC violin of n_genes / counts / pct_mt. Used to set MAD thresholds; not cell-type evidence.",
    "scatter": "QC scatter of counts vs genes or pct_mt. Inspect empty droplets, doublets, mito outliers.",
    "umap": "UMAP of the neighborhood graph. Clustering is not performed on UMAP; mixing is not integration proof.",
    "markers": "Exploratory Wilcoxon cluster markers, not a between-group result (use pseudobulk + FDR).",
    "annotation": "Annotation view. Read together with the dual-validation table.",
    "other": "Generated figure. No extra interpretation.",
}


def _fig_section(artifacts: dict, lang: str) -> str:
    caps = artifacts.get("figure_captions") or []
    if not caps:
        executed = any((p or {}).get("executed") for p in (artifacts.get("phases") or {}).values())
        if executed:
            return "- 未捕获图像文件。\n" if lang != "en" else "- No figure files captured.\n"
        return "- 未执行。图未生成。\n" if lang != "en" else "- Not executed. No figures.\n"
    lines = []
    for c in caps:
        cap = c.get("caption") or ""
        if lang == "en":
            cap = CAPTION_EN.get(c.get("kind") or "other", cap)
        lines.append(f"- `{c.get('path')}` — {cap}")
    return "\n".join(lines) + "\n"


def _metrics_table(artifacts: dict, lang: str) -> str:
    m = artifacts.get("metrics") or {}
    if not m:
        return ("未执行，无过滤数字。" if lang != "en" else "Not executed; no filter counts.")
    keys = [
        ("n_before", "过滤前细胞数" if lang != "en" else "cells before"),
        ("n_after", "过滤后细胞数" if lang != "en" else "cells after"),
        ("pct_removed", "移除比例 %" if lang != "en" else "% removed"),
        ("nmads", "MAD n"),
        ("qc_method", "QC 方法" if lang != "en" else "QC method"),
        ("imputation", "插补" if lang != "en" else "imputation"),
        ("ambient", "ambient RNA"),
        ("doublet_rate", "双细胞比例" if lang != "en" else "doublet rate"),
        ("resolution", "Leiden resolution"),
        ("n_clusters", "簇数" if lang != "en" else "clusters"),
        ("integrator", "整合方法" if lang != "en" else "integrator"),
        ("ilisi", "iLISI"),
        ("kbet", "kBET"),
        ("pca_batch_r2", "PCA 批次 R²"),
        ("batch_cluster_dominance", "cluster 内主导批次比例"),
        ("celltypist_model", "CellTypist 模型"),
        ("seed", "seed"),
    ]
    rows = ["| 指标 | 值 |", "|---|---|"] if lang != "en" else ["| metric | value |", "|---|---|"]
    for k, label in keys:
        if k in m and m[k] is not None:
            rows.append(f"| {label} | {m[k]} |")
    return "\n".join(rows) if len(rows) > 2 else ("无结构化指标。" if lang != "en" else "No structured metrics.")


def _review_block(review: dict | None, title: str) -> str:
    review = review or {}
    issues = review.get("issues") or []
    issue_txt = "；".join(issues) if issues else "无"
    return f"### {title}\n\n- 通过: {review.get('passed')}\n- 问题: {issue_txt}\n"


def _one_lang(state: dict, lang: str) -> str:
    meta = state.get("metadata") or {}
    plan = state.get("plan") or {}
    qc = state.get("qc_strategy") or {}
    ann = state.get("annotation_plan") or {}
    artifacts = state.get("artifacts") or {}
    zh = lang != "en"
    title = "scAgent 单细胞分析报告" if zh else "scAgent scRNA-seq report"
    lines = [
        f"# {title}",
        "",
        f"- {'日期' if zh else 'date'}: {date.today().isoformat()}",
        f"- {'任务' if zh else 'query'}: {state.get('user_query') or ''}",
        f"- species: {meta.get('species')} | platform: {meta.get('platform')} | tissue: {meta.get('tissue')} | n_samples: {meta.get('n_samples')}",
        f"- status: {state.get('status')}",
        f"- skills fingerprint: {artifacts.get('skills_fingerprint')}",
        f"- seed: {(artifacts.get('metrics') or {}).get('seed', 0)}",
        "",
        f"## {'数据概况' if zh else 'Dataset'}",
        "",
        f"- n_cells: {meta.get('n_cells')} | n_genes: {meta.get('n_genes')} | exists: {meta.get('exists')}",
        f"- notes: {'; '.join(meta.get('notes') or []) or ('无' if zh else 'none')}",
        "",
        f"## {'分析路线' if zh else 'Route'}",
        "",
        plan.get("narrative") or "",
        "",
        f"- skills: {', '.join(plan.get('skills') or [])}",
        f"- route: {' → '.join(plan.get('route') or [])}",
        f"- integrator: {plan.get('integrator')}",
        f"- imputation: {plan.get('imputation') or qc.get('imputation')}",
        f"- qc method: {qc.get('method')}",
        f"- mito note: {qc.get('pct_mt_note') or ('无' if zh else 'none')}",
        f"- risks: {'; '.join(plan.get('risks') or []) or ('无' if zh else 'none')}",
        "",
    ]
    if plan.get("r_degraded"):
        lines += [
            "## R 降级" if zh else "## R degraded",
            "",
            "未生成 Seurat 代码。现有 SOP 为 Scanpy。改用 `--language python` 或按规划在 R 中手工实现。",
            "",
        ]
    lines += [
        f"## {'QC 决策与数字' if zh else 'QC decisions'}",
        "",
        qc.get("protocol") or "",
        "",
        _metrics_table(artifacts, lang),
        "",
        f"## {'整合理由' if zh else 'Integration'}",
        "",
        f"{'方法' if zh else 'method'}: {plan.get('integrator') or ('未整合 / 单样本' if zh else 'none / single sample')}",
        f"{'可选模块' if zh else 'modules'}: Harmony / scVI / Scanorama(CCA-like)；Seurat CCA 需 R，不自动生成。",
        "",
        f"## {'聚类参数' if zh else 'Clustering'}",
        "",
        f"- Leiden resolution: {(artifacts.get('metrics') or {}).get('resolution') or plan.get('resolution') or 'adaptive/silhouette or 0.6'}",
        f"- n_clusters: {(artifacts.get('metrics') or {}).get('n_clusters', '未执行' if zh else 'not executed')}",
        "",
        f"## {'注释证据' if zh else 'Annotation evidence'}",
        "",
        ann.get("instructions") or "",
        "",
        f"- dual validation: {ann.get('dual_validation')}",
        f"- hierarchical: {ann.get('hierarchical')}",
        f"- forbid single gene: {ann.get('forbid_single_gene')}",
        f"- catalog: {ann.get('catalog_tissue')} ({ann.get('n_cell_types')} types)",
        f"- {'层级' if zh else 'lineage'}: Immune → T cell → CD8 T → Tex（见 cell_type_l1..）",
        "",
        f"## {'审查结论' if zh else 'Review'}",
        "",
        format_review_card(state.get("review_publication") or publication_review(state), lang),
        _review_block(state.get("review_qc"), "QC 代码/执行" if zh else "QC code/execution"),
        _review_block(state.get("review_downstream"), "聚类/注释 代码/执行" if zh else "Downstream code/execution"),
        "",
        f"## {'图' if zh else 'Figures'}",
        "",
        _fig_section(artifacts, lang),
        "",
        f"## {'警告' if zh else 'Warnings'}",
        "",
    ]
    warns = artifacts.get("warnings") or []
    lines.append("\n".join(f"- {w}" for w in warns) if warns else ("- 无" if zh else "- none"))
    if zh:
        lines += [
            "",
            "## 局限与下一步",
            "",
            "- 统计推断以生物学重复为单位；Wilcoxon 仅为探索。",
            "- UMAP 是可视化。注释是分层证据。",
            "- 未执行的步骤不得写成结果。",
            "- HITL: `--interrupt` 后检查 QC，再用 `--annotate-only --resume` 继续。",
            "",
            "## 可复现",
            "",
        ]
    else:
        lines += [
            "",
            "## Limitations",
            "",
            "- Infer at the biological-replicate level; Wilcoxon is exploratory.",
            "- UMAP is visualization. Annotation is layered evidence.",
            "- Do not describe plots that were not generated.",
            "- HITL: `--interrupt` then `--annotate-only --resume`.",
            "",
            "## Reproducibility",
            "",
        ]
    lines += [
        f"- thread_id: {state.get('thread_id')}",
        f"- jail: {(state.get('execution_qc') or {}).get('jail')}",
        f"- snapshots: {', '.join((state.get('execution_qc') or {}).get('snapshots') or []) or ('无' if zh else 'none')}",
        "- structured log: outputs/run_log.json",
        "- provenance: outputs/memory.yaml（步骤与参数，不是聊天）",
        "",
    ]
    mem = state.get("analysis_memory")
    if mem:
        from agents.memory import dump_memory_yaml

        lines += ["```yaml", dump_memory_yaml(mem).rstrip(), "```", ""]
    return "\n".join(lines)


def render_report(state: dict) -> str:
    lang = state.get("report_lang") or "zh"
    if lang == "both":
        return _one_lang(state, "zh") + "\n\n---\n\n" + _one_lang(state, "en")
    return _one_lang(state, "en" if lang == "en" else "zh")


def render_html(state: dict) -> str:
    md = render_report(state)
    arts = state.get("artifacts") or {}
    imgs = []
    for cap in arts.get("figure_captions") or []:
        p = cap.get("path") or ""
        if p:
            kind = cap.get("kind") or "figure"
            cap_txt = cap.get("caption") or p
            imgs.append(
                f'<figure><img src="{p}" alt="{kind}" style="max-width:100%"/>'
                f"<figcaption>{cap_txt}</figcaption></figure>"
            )
    body = md.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'><title>scAgent report</title>"
        "<style>body{font-family:sans-serif;max-width:920px;margin:2rem auto}"
        "pre{white-space:pre-wrap}img{max-width:100%;height:auto}figure{margin:1.5rem 0}</style>"
        f"</head><body><pre>{body}</pre>{''.join(imgs)}</body></html>"
    )


def write_run_log(state: dict, out_dir) -> None:
    import json
    from pathlib import Path

    from agents.memory import build_memory
    from scagent.config import analysis_params

    payload = {
        "thread_id": state.get("thread_id"),
        "status": state.get("status"),
        "query": state.get("user_query"),
        "params": analysis_params(),
        "provenance": state.get("analysis_memory") or build_memory(state),
        "skills": state.get("skills_used") or (state.get("plan") or {}).get("skills"),
        "route": (state.get("plan") or {}).get("route"),
        "integrator": (state.get("plan") or {}).get("integrator"),
        "metrics": (state.get("artifacts") or {}).get("metrics"),
        "review_qc": {
            "passed": (state.get("review_qc") or {}).get("passed"),
            "issues": (state.get("review_qc") or {}).get("issues"),
            "issue_records": (state.get("review_qc") or {}).get("issue_records"),
        },
        "review_downstream": {
            "passed": (state.get("review_downstream") or {}).get("passed"),
            "issues": (state.get("review_downstream") or {}).get("issues"),
            "issue_records": (state.get("review_downstream") or {}).get("issue_records"),
        },
        "jail": (state.get("execution_qc") or {}).get("jail"),
        "snapshots": (state.get("execution_qc") or {}).get("snapshots"),
    }
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    (Path(out_dir) / "run_log.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
