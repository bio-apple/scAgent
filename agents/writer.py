from __future__ import annotations

from datetime import date

from agents.reviewer import format_review_card, publication_review


CAPTION_EN = {
    "violin": "QC violin of n_genes / counts / pct_mt. Used to set MAD thresholds; not cell-type evidence.",
    "scatter": "QC scatter of counts vs genes or pct_mt. Inspect empty droplets, doublets, mito outliers.",
    "umap": "UMAP of the neighborhood graph. Clustering is not performed on UMAP; mixing is not integration proof.",
    "batch_pca_before": "PCA before integration, colored by batch. Mixing is not proof of success.",
    "batch_pca_after": "PCA/latent after integration, colored by batch. Read with iLISI/kBET.",
    "batch_umap_before": "UMAP from uncorrected PCA, colored by batch. Diagnostic only.",
    "batch_umap_after": "UMAP after integration, colored by batch. Mixing is not integration proof.",
    "markers": "Exploratory cluster markers (Wilcoxon/t-test/MAST). Not a between-group result (use pseudobulk + FDR).",
    "marker_heatmap": "Exploratory cluster marker heatmap (top genes × clusters). Not a group-level DE result.",
    "volcano": "Pseudobulk group DE volcano (logFC vs -log10 p). Interpret with FDR.",
    "pathway_bubble": "Pathway enrichment bubble plot (ORA/GSEA; -log10 p and overlap size).",
    "annotation": "Annotation view. Read together with the dual-validation table.",
    "paga": "PAGA cluster graph. Connectivity is not proof of a biological fate.",
    "pseudotime": "Diffusion pseudotime on the neighborhood graph. Exploratory axis, not a clock.",
    "gene_trends": "Gene expression vs pseudotime. Dynamic trends, not a mechanism.",
    "velocity": "RNA velocity embedding. Requires spliced/unspliced; check phase portraits.",
    "other": "Generated figure. No extra interpretation.",
}


def _fig_section(artifacts: dict, lang: str) -> str:
    all_caps = artifacts.get("figure_captions") or []
    caps = [c for c in all_caps if not str(c.get("kind") or "").startswith("batch_")]
    if not all_caps:
        executed = any((p or {}).get("executed") for p in (artifacts.get("phases") or {}).values())
        if executed:
            return "- 未捕获图像文件。\n" if lang != "en" else "- No figure files captured.\n"
        return "- 未执行。图未生成。\n" if lang != "en" else "- Not executed. No figures.\n"
    if not caps:
        return (
            "- 批次诊断图见上方 Integration。\n"
            if lang != "en"
            else "- Batch diagnostics are in Integration above.\n"
        )
    lines = []
    for c in caps:
        cap = c.get("caption") or ""
        if lang == "en":
            cap = CAPTION_EN.get(c.get("kind") or "other", cap)
        path = c.get("path") or ""
        lines.append(f"![{cap}]({path})")
        lines.append(f"- `{path}` — {cap}")
    return "\n".join(lines) + "\n"


def _batch_fig_section(artifacts: dict, lang: str) -> str:
    caps = [c for c in (artifacts.get("figure_captions") or []) if str(c.get("kind") or "").startswith("batch_")]
    if not caps:
        return ""
    zh = lang != "en"
    lines = [f"### {'校正前后批次着色' if zh else 'Batch diagnostics (before / after)'}", ""]
    order = ("batch_pca_before", "batch_pca_after", "batch_umap_before", "batch_umap_after")
    ranked = sorted(caps, key=lambda c: order.index(c.get("kind")) if c.get("kind") in order else 99)
    for c in ranked:
        cap = c.get("caption") or ""
        if not zh:
            cap = CAPTION_EN.get(c.get("kind") or "other", cap)
        path = c.get("path") or ""
        lines.append(f"![{cap}]({path})")
        lines.append(f"*{cap}*")
        lines.append("")
    return "\n".join(lines)


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
        ("doublet_rate_high_conf", "高置信双细胞比例" if lang != "en" else "high-conf doublet rate"),
        ("doublet_rate_low_conf", "低置信双细胞比例" if lang != "en" else "low-conf doublet rate"),
        ("doublet_n_high_conf", "高置信双细胞数" if lang != "en" else "high-conf doublets"),
        ("doublet_n_low_conf", "低置信双细胞数" if lang != "en" else "low-conf doublets"),
        ("doublet_filter", "双细胞过滤模式" if lang != "en" else "doublet filter"),
        ("doublet_methods", "双细胞方法" if lang != "en" else "doublet methods"),
        ("deg_engine", "DEG 后端" if lang != "en" else "DEG engine"),
        ("deg_n_sig", "DEG 显著基因数" if lang != "en" else "DEG n_sig"),
        ("deg_n_overlap", "DEG 交叉验证重叠" if lang != "en" else "DEG overlap"),
        ("marker_method", "cluster marker 方法" if lang != "en" else "marker method"),
        ("marker_n_overlap", "marker 交叉验证重叠" if lang != "en" else "marker overlap"),
        ("resolution", "Leiden resolution"),
        ("n_clusters", "簇数" if lang != "en" else "clusters"),
        ("integrator", "整合方法" if lang != "en" else "integrator"),
        ("ilisi", "iLISI"),
        ("kbet", "kBET"),
        ("pca_batch_r2", "PCA 批次 R²"),
        ("batch_cluster_dominance", "cluster 内主导批次比例"),
        ("trajectory_verdict", "轨迹判断" if lang != "en" else "trajectory verdict"),
        ("trajectory_methods", "轨迹方法" if lang != "en" else "trajectory methods"),
        ("celltypist_model", "CellTypist 模型"),
        ("seed", "seed"),
    ]
    rows = ["| 指标 | 值 |", "|---|---|"] if lang != "en" else ["| metric | value |", "|---|---|"]
    for k, label in keys:
        if k in m and m[k] is not None:
            rows.append(f"| {label} | {m[k]} |")
    return "\n".join(rows) if len(rows) > 2 else ("无结构化指标。" if lang != "en" else "No structured metrics.")


def _doublet_tier_note(artifacts: dict, lang: str) -> str:
    m = artifacts.get("metrics") or {}
    if m.get("doublet_n_high_conf") is None and m.get("doublet_rate_high_conf") is None:
        return ""
    zh = lang != "en"
    filt = m.get("doublet_filter") or "high_conf"
    lines = [
        "",
        f"**{'双细胞置信度分级' if zh else 'Doublet confidence tiers'}**",
        "",
        f"- {'高置信' if zh else 'high_conf'} (`doublet_high_conf`): {m.get('doublet_n_high_conf', '—')} ({m.get('doublet_rate_high_conf', '—')})",
        f"- {'低置信' if zh else 'low_conf'} (`doublet_low_conf`): {m.get('doublet_n_low_conf', '—')} ({m.get('doublet_rate_low_conf', '—')})",
        f"- {'过滤模式' if zh else 'filter mode'}: `{filt}`"
        + (
            "（保守：仅移除高置信；低置信保留供人工判断）"
            if zh and filt == "high_conf"
            else (" (conservative: remove high_conf only)" if not zh and filt == "high_conf" else "")
        )
        + (
            "（严格：高+低均移除）" if zh and filt == "all" else (" (strict: remove high+low)" if not zh and filt == "all" else "")
        ),
        "",
    ]
    return "\n".join(lines)


def _manifest_provenance_section(state: dict, lang: str) -> str:
    """Summarize workspace/run_manifest.json for the publication report."""
    import json

    from scagent.config import load_config, resolve_path

    path = resolve_path(load_config(), "workspace") / "run_manifest.json"
    if not path.is_file():
        return ""
    try:
        m = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return ""
    zh = lang != "en"
    env = m.get("environment") or {}
    seed_p = m.get("seed_propagation") or {}
    lines = [
        f"### {'运行清单 (run_manifest)' if zh else 'Run manifest'}",
        "",
        f"- `environment.hash`: `{env.get('hash', '—')}` ({', '.join(env.get('sources') or []) or 'n/a'})",
        f"- `seed`: {m.get('seed', '—')} | master: {seed_p.get('master_seed', '—')}",
    ]
    steps = seed_p.get("steps") or {}
    if steps:
        core = ("hvg", "pca", "neighbors", "leiden", "umap")
        lines.append(
            "- "
            + ("随机步骤" if zh else "stochastic steps")
            + ": "
            + ", ".join(f"{k}={steps.get(k)}" for k in core if k in steps)
        )
    prov = m.get("step_provenance") or []
    if prov:
        lines.append("")
        lines.append(f"**{'步骤 I/O' if zh else 'Step I/O'}**")
        lines.append("")
        lines.append("| step | input n_obs | output n_obs | obs cols (out) |")
        lines.append("|---|---:|---:|---|")
        for row in prov:
            inp = row.get("input") or {}
            out = row.get("output") or {}
            obs_out = out.get("obs_columns") or []
            obs_txt = ", ".join(obs_out[:8]) + ("…" if len(obs_out) > 8 else "")
            lines.append(
                f"| {row.get('step', '—')} | {inp.get('n_obs', '—')} | {out.get('n_obs', '—')} | {obs_txt or '—'} |"
            )
    lines.append("")
    return "\n".join(lines)


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
        f"- scAgent: {artifacts.get('scagent_version') or ''}",
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
        f"- best_practices: {', '.join(plan.get('best_practices') or [])}",
        f"- route: {' → '.join(plan.get('route') or [])}",
        f"- integrator: {plan.get('integrator')}",
        f"- imputation: {plan.get('imputation') or qc.get('imputation')}",
        f"- qc method: {qc.get('method')}",
        f"- mito note: {qc.get('pct_mt_note') or ('无' if zh else 'none')}",
        f"- risks: {'; '.join(plan.get('risks') or []) or ('无' if zh else 'none')}",
        "",
    ]
    from agents.literature import format_literature_report_block

    lines += [
        format_literature_report_block(
            plan=plan,
            qc=qc,
            ann=ann,
            interpret=state.get("interpretation_plan") or {},
            lang=lang,
        ).rstrip(),
        "",
    ]
    tr = plan.get("tool_route")
    if tr:
        from scagent.tool_router import format_route_table

        lines += [format_route_table(tr, lang=lang), ""]
    agents = plan.get("agents") or []
    if agents:
        lines += [
            f"## {'多智能体分工' if zh else 'Agents'}",
            "",
            f"{'协作' if zh else 'collaboration'}: {plan.get('collaboration') or 'multi-agent'}",
        ]
        for a in agents:
            lines.append(f"- **{a.get('name')}**: {a.get('charge')}")
        lines.append("")
    if plan.get("r_degraded"):
        lines += [
            "## R 降级" if zh else "## R degraded",
            "",
            "已写出双重格式 `outputs/analysis.Rmd`（Seurat 可运行块）。scAgent 不执行 R kernel。现有 SOP 为 Scanpy。",
            "",
        ]
    lines += [
        f"## {'QC 决策与数字' if zh else 'QC decisions'}",
        "",
        qc.get("protocol") or "",
        "",
        _metrics_table(artifacts, lang),
        _doublet_tier_note(artifacts, lang),
        "",
        f"## {'整合理由' if zh else 'Integration'}",
        "",
        f"{'方法' if zh else 'method'}: {plan.get('integrator') or ('未整合 / 单样本' if zh else 'none / single sample')}",
        f"{'决策' if zh else 'decision'}: {plan.get('integrator_reason') or plan.get('skip_integration_reason') or ('—' if zh else 'n/a')}",
        f"{'批次列' if zh else 'batch key'}: {plan.get('sample_key') or (state.get('metadata') or {}).get('sample_key') or '—'}",
        f"{'样本数' if zh else 'n_samples'}: {plan.get('n_samples') or (state.get('metadata') or {}).get('n_samples') or '—'}",
        f"{'可选模块' if zh else 'modules'}: Harmony（auto 默认）/ scVI（≥10 万细胞或 ≥8 样本）/ Scanorama(cca) / BBKNN；Seurat CCA 需 R，不自动生成。",
        f"- iLISI: {(artifacts.get('metrics') or {}).get('ilisi', '—' if zh else 'n/a')} | kBET: {(artifacts.get('metrics') or {}).get('kbet', '—' if zh else 'n/a')} | PCA-R²: {(artifacts.get('metrics') or {}).get('pca_batch_r2', '—' if zh else 'n/a')}",
        "",
        (
            ("BBKNN 改邻居图、不产生校正后 PCA；对比看 UMAP 校正前后 + iLISI/kBET。UMAP 混匀不是整合成功的证据。" if zh else "BBKNN edits the graph, not PCA; compare UMAP before/after with iLISI/kBET. Mixing is not proof.")
            if plan.get("integrator") == "bbknn"
            else ("禁止把 UMAP 混匀当作整合成功；以 iLISI / kBET / PCA 批次 R² 与校正前后图为准。" if zh else "UMAP mixing is not integration success; use iLISI/kBET/PCA-R² and before/after plots.")
        ),
        "",
        _batch_fig_section(artifacts, lang),
        f"## {'聚类参数' if zh else 'Clustering'}",
        "",
        f"- Leiden resolution: {(artifacts.get('metrics') or {}).get('resolution') or plan.get('resolution') or 'adaptive/silhouette or 0.6'}",
        f"- n_clusters: {(artifacts.get('metrics') or {}).get('n_clusters', '未执行' if zh else 'not executed')}",
        "",
        f"## {'轨迹与细胞命运' if zh else 'Trajectory / fate'}",
        "",
        (
            f"- {'判断' if zh else 'verdict'}: {(artifacts.get('metrics') or {}).get('trajectory_verdict') or ('未执行' if zh else 'not run')}"
        ),
        f"- methods: {(artifacts.get('metrics') or {}).get('trajectory_methods') or ('—' if zh else 'n/a')}",
        f"- confidence: {(artifacts.get('metrics') or {}).get('trajectory_confidence') or ('—' if zh else 'n/a')}",
        (
            "- 离散群体不强行拟合命运轴。Palantir/DPT 是探索性分化轴；scVelo 需要 spliced/unspliced；Monocle3 需 R。"
            if zh
            else "- Do not force a fate axis on discrete clusters. Palantir/DPT are exploratory; scVelo needs spliced/unspliced; Monocle3 needs R."
        ),
        "",
    ]
    hitl_mt = state.get("hitl_mt") or {}
    hitl_res = state.get("hitl_resolution") or {}
    if hitl_mt or hitl_res:
        lines += [
            f"## {'人工确认（HITL）' if zh else 'Human-in-the-loop'}",
            "",
        ]
        if hitl_mt:
            rec = hitl_mt.get("recommended")
            chosen = state.get("qc_choice") or ("自动推荐" if zh else "auto recommended")
            lines.append(f"- {'线粒体' if zh else 'MT'}: {chosen}（{('推荐' if zh else 'recommended')}={rec}）→ `outputs/decisions/mt.html`")
            for o in hitl_mt.get("options") or []:
                lines.append(f"  - `{o.get('id')}` MAD n={o.get('nmads')}: {o.get('reason')}")
        if hitl_res:
            chosen = state.get("resolution_choice") or state.get("resolution") or ("待确认" if zh else "pending")
            lines.append(f"- Leiden: {chosen} → `outputs/decisions/resolution.html`")
            for o in hitl_res.get("options") or []:
                lines.append(f"  - `{o.get('id')}` resolution={o.get('resolution')}: {o.get('reason')}")
        lines.append("")
    lines += [
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
        f"## {'生物学解读' if zh else 'Biological interpretation'}",
        "",
        (state.get("interpretation_plan") or {}).get("instructions") or ("未运行 Interpretation Agent。" if zh else "Interpretation agent did not run."),
        "",
        f"- method: {(state.get('interpretation_plan') or {}).get('method') or ('—' if zh else 'n/a')}",
        f"- gene sets: {(state.get('interpretation_plan') or {}).get('gene_sets') or ('—' if zh else 'n/a')}",
        f"- GSVA: {(state.get('interpretation_plan') or {}).get('gsva_note') or ('未执行' if zh else 'not run')}",
        f"- enrichment engine: {(artifacts.get('metrics') or {}).get('enrichment_engine', '—' if zh else 'n/a')}",
        f"- n_pathway_terms: {(artifacts.get('metrics') or {}).get('n_pathway_terms', '—' if zh else 'n/a')}",
        "",
        f"## {'证据链（细胞状态断言）' if zh else 'Evidence chain (cell-state claims)'}",
        "",
        (
            "断言「某 cluster 处于某状态」必须同时给出高表达 marker、通路 GO/id 的 p 值、以及 PubMed DOI/PMID。"
            "这是支持性证据链，不是干预因果。缺一则不得写成结果。"
            if zh
            else "A claim such as 'this cluster is exhausted T cells' requires markers, a pathway p-value, and a PubMed DOI/PMID. Supporting chain, not interventional causality."
        ),
        "",
    ]
    from scagent.evidence import render_evidence_markdown

    lines.append(render_evidence_markdown(artifacts.get("evidence_chains"), lang=lang).rstrip())
    lines += [
        "",
        f"## {'审查结论' if zh else 'Review'}",
        "",
        format_review_card(state.get("review_publication") or publication_review(state), lang),
        _review_block(state.get("review_qc"), "QC 代码/执行" if zh else "QC code/execution"),
        _review_block(state.get("review_downstream"), "聚类/DEG 代码/执行" if zh else "Cluster/DEG code/execution"),
        "",
    ]
    from scagent.publication_figures import render_publication_figure_inventory_markdown

    lines += [
        f"## {'发表级图表清单' if zh else 'Publication figure checklist'}",
        "",
        render_publication_figure_inventory_markdown(state, lang=lang).rstrip(),
        "",
        f"## {'图' if zh else 'Figures'}",
        "",
        _fig_section(artifacts, lang),
        "",
        f"## {'代码-结果（双重输出）' if zh else 'Code–result dual output'}",
        "",
    ]
    from scagent.dual import render_dual_markdown

    lines.append(render_dual_markdown(state, lang=lang, heading=False).rstrip())
    lines += [
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
            "- HITL: `--interrupt` 后打开 `outputs/decisions/*.html`，用 `scagent confirm mt|resolution <选项>` 确认再继续。",
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
            "- HITL: `--interrupt`, inspect `outputs/decisions/*.html`, then `scagent confirm mt|resolution <choice>`.",
            "",
            "## Reproducibility",
            "",
        ]
    manifest_blk = _manifest_provenance_section(state, lang)
    if manifest_blk:
        lines.append(manifest_blk.rstrip())
    lines += [
        f"- thread_id: {state.get('thread_id')}",
        f"- jail: {(state.get('execution_qc') or {}).get('jail')}",
        f"- snapshots: {', '.join((state.get('execution_qc') or {}).get('snapshots') or []) or ('无' if zh else 'none')}",
        "- structured log: outputs/run_log.json",
        "- manifest: workspace/run_manifest.json",
        f"- provenance: outputs/memory.yaml ({'步骤与参数' if zh else 'steps + params'})",
        f"- dual: outputs/dual.md ({'每阶段 [结论]+[代码]' if zh else 'conclusion + code per phase'})",
        "- notebook: outputs/analysis.ipynb",
        "- viewer: outputs/viewer.html",
        "- decisions: outputs/decisions/mt.html, resolution.html",
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


def stage_report_figures(artifacts: dict, out_dir) -> dict:
    """Copy figure files next to the report so markdown/HTML `figures/*.png` paths resolve."""
    import shutil
    from pathlib import Path

    out_dir = Path(out_dir)
    dest = out_dir / "figures"
    dest.mkdir(parents=True, exist_ok=True)
    arts = dict(artifacts or {})
    new_caps: list[dict] = []
    seen: set[str] = set()
    staged: list[str] = []
    for cap in arts.get("figure_captions") or []:
        item = dict(cap)
        src = Path(str(item.get("path") or ""))
        if src.is_file():
            name = src.name
            target = dest / name
            if src.resolve() != target.resolve():
                shutil.copy2(src, target)
            item["path"] = f"figures/{name}"
        rel = str(item.get("path") or "")
        if not rel or rel in seen:
            continue
        seen.add(rel)
        new_caps.append(item)
        staged.append(rel)
    arts["figure_captions"] = new_caps
    if staged:
        arts["figures"] = list(dict.fromkeys([*staged, *(arts.get("figures") or [])]))
    return arts


def render_html(state: dict) -> str:
    md = render_report(state)
    arts = state.get("artifacts") or {}
    integ, other = [], []
    for cap in arts.get("figure_captions") or []:
        p = cap.get("path") or ""
        if not p:
            continue
        kind = cap.get("kind") or "figure"
        cap_txt = cap.get("caption") or p
        tag = (
            f'<figure><img src="{p}" alt="{kind}" style="max-width:100%"/>'
            f"<figcaption>{cap_txt}</figcaption></figure>"
        )
        (integ if str(kind).startswith("batch_") else other).append(tag)
    gallery = ""
    if integ:
        gallery += "<h2>整合诊断（校正前后批次着色）</h2><div class='integ-gallery'>" + "".join(integ) + "</div>"
    if other:
        gallery += "<h2>Figures</h2>" + "".join(other)
    body = md.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return (
        "<!DOCTYPE html><html><head><meta charset='utf-8'><title>scAgent report</title>"
        "<style>body{font-family:sans-serif;max-width:920px;margin:2rem auto}"
        "pre{white-space:pre-wrap}img{max-width:100%;height:auto}figure{margin:1.5rem 0}"
        ".integ-gallery figure{display:inline-block;width:48%;vertical-align:top}</style>"
        f"</head><body><pre>{body}</pre>{gallery}</body></html>"
    )


def write_run_log(state: dict, out_dir) -> None:
    import json
    from pathlib import Path

    from agents.memory import build_memory
    from scagent.config import analysis_params
    from scagent.export_nb import package_versions
    from scagent.publication_figures import build_publication_figure_inventory

    payload = {
        "thread_id": state.get("thread_id"),
        "status": state.get("status"),
        "query": state.get("user_query"),
        "params": analysis_params(),
        "versions": package_versions(),
        "provenance": state.get("analysis_memory") or build_memory(state),
        "skills": state.get("skills_used") or (state.get("plan") or {}).get("skills"),
        "route": (state.get("plan") or {}).get("route"),
        "integrator": (state.get("plan") or {}).get("integrator"),
        "integrator_reason": (state.get("plan") or {}).get("integrator_reason"),
        "metrics": (state.get("artifacts") or {}).get("metrics"),
        "publication_figures": build_publication_figure_inventory(state),
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
        "executor": (state.get("execution_qc") or {}).get("jail"),
        "jail": (state.get("execution_qc") or {}).get("jail"),
        "snapshots": (state.get("execution_qc") or {}).get("snapshots"),
        "snapshot_manifests": (state.get("execution_qc") or {}).get("snapshot_manifests"),
    }
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    (Path(out_dir) / "run_log.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
