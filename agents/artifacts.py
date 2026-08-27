from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

from scagent.skills_loader import list_skills

METRICS_RE = re.compile(r"^SCAGENT_METRICS:(.*)$", re.M)
WARN_RE = re.compile(r"^SCAGENT_WARN:(.*)$", re.M)

FIGURE_CAPTIONS = {
    "batch_pca_before": "校正前 PCA，按批次着色。用于看整合前批次是否分离；混匀不是成功证据。",
    "batch_pca_after": "校正后 PCA/潜空间，按批次着色。与校正前对照阅读，并同时看 iLISI/kBET。",
    "batch_umap_before": "校正前 UMAP（由未校正 PCA 计算）。仅作批次诊断，不是聚类坐标。",
    "batch_umap_after": "校正后 UMAP，按批次着色。混匀不能单独当作整合成功。",
    "violin": "QC violin：n_genes / total_counts / pct_mt 分布。用于定 MAD 阈值，不能单独当细胞类型证据。",
    "scatter": "QC scatter：counts vs genes 或 counts vs pct_mt。用于看空液滴、双细胞与线粒体离群。",
    "umap": "UMAP 可视化邻域图。聚类不在 UMAP 坐标上做；混匀不是整合成功的证据。",
    "markers": "探索性 cluster marker（Wilcoxon / t-test / MAST）。不是组间结论；组间比较需 pseudobulk + FDR。",
    "annotation": "注释验证图（参考标签或 marker）。需与 dual validation 表一起读。",
    "paga": "PAGA 簇图。连接强度不是生物学命运的证明。",
    "pseudotime": "扩散伪时间。探索性分化轴，不是时钟。",
    "gene_trends": "基因随伪时间的趋势。动态相关，不是机制结论。",
    "velocity": "RNA velocity 嵌入。需要 spliced/unspliced；并检查 phase portrait。",
}


def skills_fingerprint() -> str:
    h = hashlib.sha256()
    for skill in list_skills():
        h.update(skill.name.encode())
        h.update(skill.path.read_bytes())
    return h.hexdigest()[:12]


def parse_metrics(stdout: str, stderr: str) -> tuple[dict[str, Any], list[str]]:
    metrics: dict[str, Any] = {}
    m = METRICS_RE.search(stdout or "")
    if m:
        try:
            metrics = json.loads(m.group(1))
        except json.JSONDecodeError:
            metrics = {"raw": m.group(1)[:500]}
    warns = [w.strip() for w in WARN_RE.findall(stdout or "")]
    warns += [w.strip() for w in WARN_RE.findall(stderr or "")]
    return metrics, warns


def collect_workspace(workspace: Path, phase: str, execution: dict[str, Any] | None = None) -> dict[str, Any]:
    execution = execution or {}
    fig_dir = workspace / "figures"
    figures = sorted(str(p) for p in fig_dir.glob("*") if p.is_file()) if fig_dir.exists() else []
    stdout = execution.get("stdout") or ""
    stderr = execution.get("stderr") or ""
    metrics, warns = parse_metrics(stdout, stderr)
    h5ads = {
        "qc": str(workspace / "adata_qc.h5ad") if (workspace / "adata_qc.h5ad").exists() else None,
        "processed": str(workspace / "adata_processed.h5ad")
        if (workspace / "adata_processed.h5ad").exists()
        else None,
    }
    captions = []
    for path in figures:
        low = Path(path).name.lower()
        kind = "other"
        for key in sorted(FIGURE_CAPTIONS, key=len, reverse=True):
            if key in low:
                kind = key
                break
        captions.append(
            {
                "path": path,
                "kind": kind,
                "caption": FIGURE_CAPTIONS.get(kind, "已生成图像。未标注类型，不解释未观测现象。"),
            }
        )
    evidence_chains = None
    chains_path = workspace / "evidence_chains.json"
    if chains_path.is_file():
        try:
            evidence_chains = json.loads(chains_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            evidence_chains = None
    out = {
        "phase": phase,
        "executed": bool(execution.get("executed")),
        "ok": bool(execution.get("ok", True)),
        "returncode": 0 if execution.get("ok", True) else 1,
        "script": execution.get("script"),
        "figures": figures,
        "figure_captions": captions,
        "h5ads": h5ads,
        "metrics": metrics,
        "warnings": warns,
        "stderr_tail": (stderr or "")[-1500:],
        "stdout_tail": (stdout or "")[-1500:],
        "evidence_chains": evidence_chains,
    }
    return out


def merge_artifacts(prev: dict[str, Any] | None, new: dict[str, Any]) -> dict[str, Any]:
    prev = dict(prev or {})
    phases = dict(prev.get("phases") or {})
    phases[new["phase"]] = new
    figures = list(dict.fromkeys([*(prev.get("figures") or []), *(new.get("figures") or [])]))
    captions = [*(prev.get("figure_captions") or []), *(new.get("figure_captions") or [])]
    warnings = [*(prev.get("warnings") or []), *(new.get("warnings") or [])]
    metrics = dict(prev.get("metrics") or {})
    metrics.update(new.get("metrics") or {})
    h5ads = dict(prev.get("h5ads") or {})
    h5ads.update(new.get("h5ads") or {})
    chains = new.get("evidence_chains") or prev.get("evidence_chains")
    return {
        "phases": phases,
        "figures": figures,
        "figure_captions": captions,
        "warnings": warnings,
        "metrics": metrics,
        "h5ads": h5ads,
        "evidence_chains": chains,
        "skills_fingerprint": prev.get("skills_fingerprint") or skills_fingerprint(),
        "scagent_version": prev.get("scagent_version") or new.get("scagent_version"),
    }
