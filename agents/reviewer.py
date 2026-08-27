from __future__ import annotations

import re
from pathlib import Path

from agents.code_schema import validate_script
from agents.common import read_prompt, run_specialist
from agents.markers import choose_celltypist_model
from agents.templates import LOCKED_END, LOCKED_START
from scagent.config import load_config
from scagent.doublets import COMPLEX_TISSUES


def _records_and_messages(records: list[dict]) -> tuple[list[dict], list[str]]:
    return records, [r["message"] for r in records]


def _add(records: list[dict], message: str, *, id: str, severity: str = "fail", source: str = "code") -> None:
    records.append({"id": id, "severity": severity, "source": source, "message": message})


def audit_code(code: str, metadata: dict | None = None, phase: str = "qc") -> dict:
    """Deterministic statistical audit. LLM cannot override hard fails."""
    metadata = metadata or {}
    text = code or ""
    low = text.lower()
    records: list[dict] = []
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
    has_pseudobulk_impl = "pseudobulk_de" in low or "get.aggregate" in low or "sc.get.aggregate" in low
    has_ref2 = any(k in low for k in ("ref2_label", "singler", "popv", "second_reference", "ref_crossval", "deg_label", "cluster_deg"))
    has_fusion = "fuse_annotation" in low or "annotation_n_agree" in low or "annotation_fusion" in low
    has_predicted_doublet = "predicted_doublet" in low or "doublet_call" in low
    has_padj = "pvals_adj" in low or "padj" in low or "fdr" in low
    tissue = str((metadata or {}).get("tissue") or "").lower()
    ct_model = choose_celltypist_model(tissue, (metadata or {}).get("species"))

    if phase == "qc":
        if not has_violin:
            _add(records, "QC 缺少 Violin/VlnPlot", id="qc.violin")
        if not has_scatter:
            _add(records, "QC 缺少 Scatter", id="qc.scatter")
        if not has_mad:
            _add(records, "QC 未使用 MAD 自适应阈值", id="qc.mad")
        if not has_locked:
            _add(records, "QC 缺少不可删除的 LOCKED QC 代码块", id="qc.locked")
        if has_mad and not mt_one_sided:
            _add(records, "线粒体 MAD 应为单侧（只滤高 pct_mt）", id="qc.mt_side")
        if not log1p_qc:
            _add(records, "calculate_qc_metrics 需 log1p=True，以便 MAD 使用 log1p_total_counts 列", id="qc.log1p")
        if re.search(r"pct[_]?mt\s*[<>=]+\s*(5|10)\b", low):
            _add(records, "使用了固定 pctMT=5/10 阈值，应改为看分布 + MAD/percentile", id="qc.hard_mt")
        if not has_predicted_doublet and "detect_doublets" not in low and "scrublet" not in low:
            _add(records, "QC 缺少双细胞检测（Scrublet / detect_doublets / doublet_call）", id="qc.doublet")
        if "scrublet skipped" in low and "predicted_doublet" not in low:
            _add(records, "Scrublet 被跳过且未写入 predicted_doublet", id="qc.doublet_skip")
        n_samp = int((metadata or {}).get("n_samples") or 1)
        if (metadata or {}).get("need_batch_correction"):
            n_samp = max(n_samp, 2)
        tiss = str((metadata or {}).get("tissue") or "").lower()
        if (n_samp > 1 or tiss in COMPLEX_TISSUES) and "detect_doublets" not in low and "scdblfinder" not in low:
            _add(records, "多样本/复杂组织缺少双细胞交叉验证（Scrublet + scDblFinder 或表达模拟）", id="qc.doublet_cross")

    if phase == "downstream":
        if metadata.get("need_batch_correction"):
            integrated = any(k in low for k in ("harmony", "scvi", "scanorama", "cca", "bbknn"))
            skipped = "skip integration" in low or "跳过整合" in text or "user disabled batch" in low
            if not integrated and not skipped:
                _add(records, "多样本未做整合，也未声明跳过理由", id="down.integrate")
            if integrated and re.search(r"umap.{0,40}(mix|mixed|混匀)", low) and "ilisi" not in low and "kbet" not in low:
                _add(records, "禁止把 UMAP 混匀当作整合成功；需 iLISI/kBET/PCA-R²", id="down.umap_mix")
        if "umap" in low and re.search(r"leiden\(.*umap|cluster.*umap", low):
            _add(records, "疑似在 UMAP 坐标上聚类", id="down.cluster_umap")
        if ct_model:
            if not has_celltypist:
                _add(records, "注释未调用 CellTypist（或未显式失败降级）", id="down.celltypist")
            if tissue not in {"pbmc", "blood", "immune"} and (
                "immune_all_low.pkl" in low or "immune_all_high.pkl" in low
            ):
                _add(records, "非免疫组织使用了 Immune_All CellTypist 模型", id="down.immune_all")
        else:
            if "immune_all_low.pkl" in low or "immune_all_high.pkl" in low:
                _add(records, "无匹配模型的组织不应默认 Immune_All", id="down.immune_all_default")
        if not has_ref2:
            _add(records, "缺少第二参考注释（cluster DE∩catalog / SingleR/popV，不能只靠 Azimuth）", id="down.ref2")
        if not has_fusion:
            _add(records, "注释缺少多证据融合（CellTypist + marker + DE overlap 表决，禁止只调用 Azimuth）", id="down.fusion")
        if "azimuth" in low and not has_fusion and "celltypist" not in low:
            _add(records, "注释只调用了 Azimuth，缺少独立 marker/DE 证据融合", id="down.azimuth_only")
        if not has_dual:
            _add(records, "注释缺少 dual validation（≥2 阳性 + ≥1 阴性 marker）", id="down.dual")
        if "cell_type_l1" not in low and "lineage" not in low:
            _add(records, "注释缺少层级字段（cell_type_l1 / lineage）", id="down.lineage")
        if metadata.get("needs_pseudobulk"):
            if not has_pseudobulk_impl or not (has_pseudobulk_note or "fdr" in low):
                _add(records, "组间比较必须 sample-level pseudobulk + FDR，不能只用 cell-level Wilcoxon", id="down.pseudobulk")
        elif "rank_genes_groups" in low or "rank_genes(" in low:
            if not has_pseudobulk_note:
                _add(records, "差异表达未声明探索-only / 组间须 pseudobulk+FDR", id="down.deg_note")
        if ("rank_genes_groups" in low or "rank_genes(" in low) and not has_padj:
            _add(records, "rank_genes_groups 必须使用 pvals_adj/FDR，不能只看 raw p", id="down.padj")
        if re.search(r"tl\.rank_genes_groups", low) and "scale(" in low and "use_raw" not in low:
            _add(records, "Wilcoxon 须 use_raw=True（或 .raw），不能在 scale 后的 X 上做", id="down.deg_scaled")
        if metadata.get("batch_condition_confounded") and any(k in low for k in ("harmony", "scvi", "scanorama", "bbknn")):
            skipped = "skip integration" in low or "跳过整合" in text or "collinear" in low
            if not skipped:
                _add(records, "样本与条件共线时不应整合，否则抹掉处理效应", id="down.confound")
        if "cell_type" in low and not has_dual:
            _add(records, "单基因或无双验证的细胞类型赋值", id="down.single_gene")
        if "trajectory" in (metadata.get("route") or []) and not any(
            k in low for k in ("dpt", "diffmap", "paga", "monocle", "palantir", "scvelo", "run_trajectory_phase", "velocity")
        ):
            _add(records, "路线含 trajectory 但未在 PCA/UMAP/Leiden 之后调用 DPT/PAGA/Palantir/scVelo", id="down.trajectory")
        schema = validate_script(text, phase="downstream")
        for rec in schema.get("issue_records") or []:
            records.append(rec)

    if phase == "qc":
        schema = validate_script(text, phase="qc")
        for rec in schema.get("issue_records") or []:
            records.append(rec)

    if not has_seed:
        _add(records, "未固定随机种子", id="seed", severity="warn")

    _, issues = _records_and_messages(records)

    hard = [x for x in issues if not x.startswith("未固定")]
    passed = not any(r.get("severity") == "fail" for r in records)

    return {
        "passed": passed,
        "issues": issues,
        "issue_records": records,
        "required_fixes": hard if not passed else [],
        "has_violin": has_violin,
        "has_scatter": has_scatter,
        "has_mad": has_mad,
        "has_locked": has_locked,
        "has_celltypist": has_celltypist,
        "has_dual": has_dual,
        "has_ref2": has_ref2,
        "has_fusion": has_fusion,
        "has_pseudobulk_impl": has_pseudobulk_impl,
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
    records: list[dict] = []
    if not execute_code or not execution.get("executed"):
        rec = {"id": "exec.skip", "severity": "warn", "source": "execution", "message": "未执行代码，结果指标未验证"}
        return {
            "passed": True,
            "skipped": True,
            "issues": [rec["message"]],
            "issue_records": [rec],
        }

    if not execution.get("ok"):
        _add(records, "执行失败（非零 returncode）", id="exec.nonzero", source="execution")
        err = (execution.get("stderr") or "")[-800:]
        if err:
            _add(records, f"stderr: {err[:300]}", id="exec.stderr", source="execution")

    metrics = artifacts.get("metrics") or {}
    phase_art = (artifacts.get("phases") or {}).get(phase) or {}
    metrics = {**metrics, **(phase_art.get("metrics") or {}), **(execution.get("metrics") or {})}

    if phase == "qc":
        pct = metrics.get("pct_removed")
        warn_pct = float(
            metadata.get("overfilter_warn_pct")
            if metadata.get("overfilter_warn_pct") is not None
            else (load_config().get("qc") or {}).get("overfilter_warn_pct")
            or 30
        )
        if pct is not None and float(pct) > warn_pct:
            _add(records, f"过度过滤：MAD 移除 {float(pct):.1f}% 细胞（阈值 {warn_pct:.0f}%）", id="exec.overfilter", source="execution")
        h5ad = (artifacts.get("h5ads") or {}).get("qc") or (phase_art.get("h5ads") or {}).get("qc")
        if execution.get("ok") and not h5ad:
            _add(records, "未写出 adata_qc.h5ad", id="exec.qc_h5ad", source="execution")
        figs = artifacts.get("figures") or phase_art.get("figures") or execution.get("figures") or []
        names = " ".join(Path(p).name.lower() for p in figs)
        if execution.get("ok"):
            if "violin" not in names:
                _add(records, "执行后缺少 violin 图", id="exec.violin", source="execution")
            if "scatter" not in names:
                _add(records, "执行后缺少 scatter 图", id="exec.scatter", source="execution")
        dstat = metrics.get("doublet_status")
        if dstat == "failed":
            _add(records, "双细胞检测未成功写入 predicted_doublet", id="exec.doublet", source="execution")
        rate = metrics.get("doublet_rate")
        max_rate = float((load_config().get("qc") or {}).get("doublet_rate_max") or 0.10)
        if rate is not None and float(rate) > max_rate:
            _add(records, f"双细胞比例 {float(rate):.1%} 超过 {max_rate:.0%}", id="exec.doublet_rate", source="execution")

    if phase == "downstream":
        h5ad = (artifacts.get("h5ads") or {}).get("processed")
        if execution.get("ok") and not h5ad:
            _add(records, "未写出 adata_processed.h5ad", id="exec.processed_h5ad", source="execution")
        mix = metrics.get("batch_cluster_dominance")
        if metadata.get("need_batch_correction"):
            integ = load_config().get("integration") or {}
            ilisi = metrics.get("ilisi")
            kbet = metrics.get("kbet")
            r2 = metrics.get("pca_batch_r2")
            if metrics.get("integration_passed") is False:
                _add(records, "整合质量未达标（iLISI/kBET/PCA-R²）", id="exec.integ", source="execution")
            elif ilisi is not None and float(ilisi) < float(integ.get("ilisi_min") or 0.8):
                _add(records, f"iLISI {float(ilisi):.3f} < {integ.get('ilisi_min')}", id="exec.ilisi", source="execution")
            elif kbet is not None and float(kbet) < float(integ.get("kbet_min") or 0.5):
                _add(records, f"kBET {float(kbet):.3f} < {integ.get('kbet_min')}", id="exec.kbet", source="execution")
            elif ilisi is None and kbet is None and r2 is not None and float(r2) > float(integ.get("pca_batch_r2_max") or 0.5):
                _add(records, f"PCA 批次 R² {float(r2):.3f} 过高", id="exec.pca_r2", source="execution")
            elif mix is not None and float(mix) >= 0.95:
                _add(records, f"整合质量可疑：cluster 内主导批次比例 {mix:.2f}（可能未混合）", id="exec.mix", source="execution")
            if execution.get("ok") and not _has_batch_diagnostic(artifacts, metrics):
                _add(
                    records,
                    "执行后缺少校正前后批次着色诊断图（PCA/UMAP）；不能只报 iLISI/kBET",
                    id="exec.integ_plot",
                    source="execution",
                )

    _, issues = _records_and_messages(records)
    passed = execution.get("ok") is True and not any(r.get("severity") == "fail" for r in records)
    return {"passed": passed, "skipped": False, "issues": issues, "issue_records": records, "metrics": metrics}


WEIGHTS = {
    "qc": 20,
    "batch_correction": 15,
    "doublet_detection": 10,
    "markers": 15,
    "deg": 10,
    "figures": 10,
    "annotation": 10,
    "evidence": 10,
}

ICONS = {"pass": "✅", "fail": "❌", "missing": "⚠️"}


def _has_batch_diagnostic(artifacts: dict, metrics: dict | None = None) -> bool:
    blob = " ".join(str(p).lower() for p in (artifacts.get("figures") or []))
    blob += " " + " ".join(str((c or {}).get("path") or "").lower() for c in (artifacts.get("figure_captions") or []))
    mets = metrics or artifacts.get("metrics") or {}
    blob += " " + " ".join(str(p).lower() for p in (mets.get("integration_plots") or []))
    return "batch_pca" in blob or "batch_umap" in blob


def _code(state: dict) -> str:
    return f"{state.get('code_qc') or ''}\n{state.get('code_downstream') or ''}"


def _item(key: str, status: str, detail: str) -> dict:
    w = WEIGHTS[key]
    if status == "pass":
        pts = w
    elif status == "missing":
        pts = w // 2
    else:
        pts = 0
    return {"key": key, "status": status, "detail": detail, "weight": w, "points": pts}


def publication_review(state: dict) -> dict:
    """Publication Reviewer card: PASS/FAIL/Missing checklist + 0–100 score."""
    meta = state.get("metadata") or {}
    plan = state.get("plan") or {}
    artifacts = state.get("artifacts") or {}
    rq = state.get("review_qc") or {}
    rd = state.get("review_downstream") or {}
    text = _code(state).lower()
    warns = " ".join(str(w) for w in (artifacts.get("warnings") or [])).lower()
    executed = bool(state.get("execute_code")) and any(
        (state.get(k) or {}).get("executed") for k in ("execution_qc", "execution_downstream")
    )
    figs = artifacts.get("figures") or []
    fig_names = " ".join(str(p).lower() for p in figs)

    items: list[dict] = []

    if not (state.get("code_qc") or rq):
        items.append(_item("qc", "missing", "QC 阶段未运行"))
    elif rq.get("passed") is False:
        items.append(_item("qc", "fail", "；".join(rq.get("issues") or ["QC 未通过"])))
    else:
        items.append(_item("qc", "pass", "Violin/Scatter/MAD/LOCKED 审查通过"))

    need_batch = bool(meta.get("need_batch_correction") or (int(meta.get("n_samples") or 1) > 1))
    mets = artifacts.get("metrics") or {}
    mix = mets.get("batch_cluster_dominance")
    ilisi = mets.get("ilisi")
    if not need_batch:
        items.append(_item("batch_correction", "pass", "单样本，无需批次校正"))
    elif any(k in text for k in ("harmony", "scvi", "scanorama", "bbknn")):
        if mets.get("integration_passed") is False or (ilisi is not None and float(ilisi) < 0.8):
            items.append(_item("batch_correction", "fail", f"整合度量未达标 iLISI={ilisi} mix={mix}"))
        elif mix is not None and float(mix) >= 0.95:
            items.append(_item("batch_correction", "fail", f"cluster 内主导批次 {mix:.2f}，可能未混合"))
        elif executed and not _has_batch_diagnostic(artifacts, mets):
            items.append(_item("batch_correction", "fail", "缺校正前后批次着色诊断图（PCA/UMAP）；不能只报 iLISI/kBET"))
        else:
            detail = f"整合={plan.get('integrator') or 'detected'} iLISI={ilisi}"
            if _has_batch_diagnostic(artifacts, mets):
                detail += "；已嵌入校正前后批次着色图"
            items.append(_item("batch_correction", "pass", detail))
    elif "skip integration" in text or plan.get("skip_integration_reason"):
        items.append(_item("batch_correction", "missing", plan.get("skip_integration_reason") or "已声明跳过整合"))
    else:
        items.append(_item("batch_correction", "fail", "多样本未做整合，也未声明跳过"))

    dstat = mets.get("doublet_status")
    engines = mets.get("doublet_methods") or []
    if isinstance(engines, str):
        engines = [e.strip() for e in engines.split(",") if e.strip()]
    engine_txt = "+".join(str(e) for e in engines) if engines else ""
    if "scrublet skipped" in warns or dstat == "failed":
        items.append(_item("doublet_detection", "missing" if dstat != "failed" else "fail", "双细胞检测未成功"))
    elif "predicted_doublet" in text or "doublet_call" in text or "scrublet" in text or "scdblfinder" in text or "detect_doublets" in text:
        rate = mets.get("doublet_rate")
        if rate is not None and float(rate) > 0.10:
            items.append(_item("doublet_detection", "fail", f"doublet_rate={float(rate):.1%}"))
        elif len(engines) >= 2 or ("detect_doublets" in text and ("scdblfinder" in text or "simulation" in text or "both" in text)):
            agree = mets.get("doublet_agreement")
            n_hi = mets.get("doublet_n_high_conf")
            n_lo = mets.get("doublet_n_low_conf")
            detail = "Scrublet + 第二方法 → doublet_call 三级"
            if n_hi is not None:
                detail += f" high={n_hi} low={n_lo}"
            if agree is not None:
                detail += f" agreement={agree}"
            if engine_txt:
                detail += f" ({engine_txt})"
            items.append(_item("doublet_detection", "pass", detail))
        else:
            n_hi = mets.get("doublet_n_high_conf")
            detail = "Scrublet → doublet_call 分级"
            if n_hi is not None:
                detail += f" high={n_hi} low={mets.get('doublet_n_low_conf')}"
            items.append(_item("doublet_detection", "pass", detail))
    else:
        items.append(_item("doublet_detection", "missing", "未检测到双细胞检测"))

    if rd.get("has_dual") or ("positive" in text and "negative" in text):
        items.append(_item("markers", "pass", "≥2 阳性 + ≥1 阴性 marker 双验证"))
    elif not (state.get("code_downstream") or rd):
        items.append(_item("markers", "missing", "注释阶段未运行"))
    else:
        items.append(_item("markers", "fail", "marker 双验证不完整"))

    need_pb = bool(plan.get("needs_pseudobulk") or meta.get("needs_pseudobulk"))
    eng = str(mets.get("deg_engine") or "")
    cv_n = mets.get("deg_n_overlap")
    marker_n = mets.get("marker_n_overlap")
    cv_bit = ""
    if cv_n is not None:
        cv_bit = f"；交叉验证 overlap={cv_n} jaccard={mets.get('deg_jaccard')}"
    elif marker_n is not None:
        cv_bit = f"；marker 交叉验证 overlap={marker_n} jaccard={mets.get('marker_jaccard')}"
    if need_pb:
        if "edger" in eng:
            items.append(_item("deg", "pass", f"sample-level pseudobulk + edgeR QL + FDR ({eng}){cv_bit}"))
        elif "deseq2" in eng:
            items.append(_item("deg", "pass", f"sample-level pseudobulk + DESeq2 + FDR ({eng}){cv_bit}"))
        elif "ttest" in eng:
            items.append(_item("deg", "pass", f"sample-level pseudobulk + t-test+BH ({eng}){cv_bit}"))
        elif "pseudobulk_de" in text or "get.aggregate" in text:
            items.append(_item("deg", "pass", f"sample-level pseudobulk + FDR（优先 edgeR/DESeq2 via rpy2）{cv_bit}"))
        else:
            items.append(_item("deg", "fail", "组间比较仍是 cell-level Wilcoxon/MAST"))
    elif "rank_genes_groups" in text or "rank_genes(" in text or "wilcoxon" in text or "t-test" in text:
        methods = mets.get("marker_methods") or mets.get("marker_method") or "wilcoxon"
        if "pseudobulk" in text:
            items.append(_item("deg", "pass", f"探索性 {methods}；组间须 pseudobulk+FDR{cv_bit}"))
        else:
            items.append(_item("deg", "fail", "DEG 未声明 pseudobulk/FDR 要求"))
    else:
        items.append(_item("deg", "missing", "未做差异表达（或仅 QC）"))

    if not executed:
        items.append(_item("figures", "missing", "未执行，图未生成"))
    else:
        from scagent.publication_figures import build_publication_figure_inventory

        pub = build_publication_figure_inventory(state)
        if pub["n_missing"]:
            items.append(
                _item(
                    "figures",
                    "fail",
                    "发表级清单缺图: " + ", ".join(pub["missing_ids"]),
                )
            )
        elif pub["n_required"]:
            items.append(_item("figures", "pass", f"发表级清单 {pub['n_present']}/{pub['n_required']} 已就绪；共 {len(figs)} 张图"))
        else:
            need_violin = bool(state.get("code_qc"))
            need_umap = bool(state.get("code_downstream"))
            missing_fig = []
            if need_violin and "violin" not in fig_names:
                missing_fig.append("violin")
            if need_violin and "scatter" not in fig_names:
                missing_fig.append("scatter")
            if need_umap and "umap" not in fig_names and "overview" not in fig_names:
                missing_fig.append("umap")
            if need_batch and not _has_batch_diagnostic(artifacts, mets):
                missing_fig.append("batch_pca/umap")
            if missing_fig:
                items.append(_item("figures", "fail", "缺图: " + ", ".join(missing_fig)))
            else:
                items.append(_item("figures", "pass", f"{len(figs)} 张图"))

    if rd.get("has_dual") and rd.get("has_fusion") and (
        rd.get("has_celltypist") or "celltypist" in text or rd.get("has_ref2") or "deg_label" in text
    ):
        if rd.get("passed") is False and any("注释" in x or "dual" in x.lower() or "融合" in x for x in (rd.get("issues") or [])):
            items.append(_item("annotation", "fail", "；".join(rd.get("issues") or ["注释证据不足"])))
        else:
            items.append(_item("annotation", "pass", "多证据融合：CellTypist + marker 双验证 + cluster DE"))
    elif not (state.get("code_downstream") or rd):
        items.append(_item("annotation", "missing", "注释阶段未运行"))
    else:
        items.append(_item("annotation", "fail", "缺少多证据融合（禁止只调用 Azimuth）"))

    chains = artifacts.get("evidence_chains") or {}
    claims = list(chains.get("claims") or mets.get("evidence_claims") or [])
    n_ok = int(chains.get("n_ok") if chains.get("n_ok") is not None else mets.get("n_claims_ok") or 0)
    n_claims = int(chains.get("n_claims") if chains.get("n_claims") is not None else (mets.get("n_claims") or len(claims)))
    interpret_ran = bool((state.get("execution_interpret") or {}).get("executed"))
    if not executed:
        items.append(_item("evidence", "missing", "未执行，细胞状态断言未附证据链"))
    elif n_claims:
        bad = [c for c in claims if isinstance(c, dict) and not c.get("ok")]
        if bad or n_ok < n_claims:
            items.append(_item("evidence", "fail", "有细胞状态断言但缺少 marker + 通路 p 值 + DOI/PMID"))
        else:
            items.append(_item("evidence", "pass", f"{n_ok}/{n_claims} 条断言具有 marker + GO p 值 + PubMed DOI"))
    elif interpret_ran or (state.get("code_interpret") and executed):
        items.append(_item("evidence", "pass", "无细胞状态断言（unknown/未注释），未写机制结论"))
    else:
        items.append(_item("evidence", "missing", "Interpretation 未写出证据链"))

    score = int(round(sum(i["points"] for i in items)))
    fails = [i for i in items if i["status"] == "fail"]
    if fails:
        verdict = "FAIL"
    elif any(i["status"] == "missing" for i in items):
        verdict = "PASS"
    else:
        verdict = "PASS"
    return {
        "items": items,
        "score": score,
        "max_score": 100,
        "verdict": verdict,
        "passed": not fails,
        "phase": "publication",
    }


def format_review_card(card: dict | None, lang: str = "zh") -> str:
    card = card or {}
    labels = {
        "qc": "QC",
        "batch_correction": "Batch correction",
        "doublet_detection": "Doublet detection",
        "markers": "Markers",
        "deg": "DEG",
        "figures": "Figures",
        "annotation": "Cell annotation evidence",
        "evidence": "Causal evidence chain",
    }
    status_txt = {"pass": "PASS", "fail": "FAIL", "missing": "Missing"}
    lines = ["### Reviewer 输出" if lang != "en" else "### Reviewer output", ""]
    for it in card.get("items") or []:
        icon = ICONS.get(it["status"], "•")
        st = status_txt.get(it["status"], it["status"])
        label = labels.get(it["key"], it["key"])
        extra = f" — {it['detail']}" if it.get("detail") else ""
        lines.append(f"{icon} **{label}:** {st}{extra}")
    lines += [
        "",
        f"**Overall score: {card.get('score', 0)} / {card.get('max_score', 100)}**"
        + (f"  ({card.get('verdict')})" if card.get("verdict") else ""),
        "",
    ]
    return "\n".join(lines)


def review_state(state: dict, phase: str | None = None) -> dict:
    phase = phase or state.get("phase") or "qc"
    code = state.get("code") or (state.get("code_qc") if phase == "qc" else state.get("code_downstream")) or ""
    meta = dict(state.get("metadata") or {})
    plan = state.get("plan") or {}
    if plan.get("needs_pseudobulk"):
        meta["needs_pseudobulk"] = True
    if plan.get("celltypist_model"):
        meta["celltypist_model"] = plan["celltypist_model"]
    if plan.get("route"):
        meta["route"] = plan["route"]
    qc = state.get("qc_strategy") or {}
    if qc.get("overfilter_warn_pct") is not None:
        meta["overfilter_warn_pct"] = qc["overfilter_warn_pct"]
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
    records = list(code_result.get("issue_records") or []) + list(exe_result.get("issue_records") or [])
    passed = bool(code_result.get("passed")) and bool(exe_result.get("passed"))
    result = {
        **code_result,
        "passed": passed,
        "issues": issues,
        "issue_records": records,
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
