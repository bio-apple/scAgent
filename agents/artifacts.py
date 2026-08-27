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
    "violin": "QC violin：n_genes / total_counts / pct_mt 分布。用于定 MAD 阈值，不能单独当细胞类型证据。",
    "scatter": "QC scatter：counts vs genes 或 counts vs pct_mt。用于看空液滴、双细胞与线粒体离群。",
    "umap": "UMAP 可视化邻域图。聚类不在 UMAP 坐标上做；混匀不是整合成功的证据。",
    "markers": "探索性 cluster marker（Wilcoxon）。不是组间结论；组间比较需 pseudobulk + FDR。",
    "annotation": "注释验证图（参考标签或 marker）。需与 dual validation 表一起读。",
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
        for key in FIGURE_CAPTIONS:
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
    return {
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
    }


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
    return {
        "phases": phases,
        "figures": figures,
        "figure_captions": captions,
        "warnings": warnings,
        "metrics": metrics,
        "h5ads": h5ads,
        "skills_fingerprint": prev.get("skills_fingerprint") or skills_fingerprint(),
    }
