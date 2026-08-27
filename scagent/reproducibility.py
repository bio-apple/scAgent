"""Environment hashing, step I/O provenance, and seed propagation metadata for run_manifest.json."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from scagent.config import analysis_params
from scagent.logutil import get_logger

log = get_logger("reproducibility")

_SEED_STEPS = ("numpy", "scanpy_settings", "hvg", "pca", "neighbors", "leiden", "umap", "scrublet", "scvi")


def seed_propagation_record(seed: int | None = None) -> dict[str, Any]:
    """Document which stochastic steps use the analysis seed."""
    s = int(seed if seed is not None else analysis_params()["seed"])
    return {
        "master_seed": s,
        "steps": {k: s for k in _SEED_STEPS},
        "pythonhashseed": "set via kernel/subprocess env at execute time",
        "note": "Set PYTHONHASHSEED before Python starts for full hash reproducibility.",
    }


def apply_analysis_seed(seed: int | None = None) -> int:
    """Apply numpy + Scanpy global seeds (call at script start)."""
    import numpy as np

    s = int(seed if seed is not None else analysis_params()["seed"])
    np.random.seed(s)
    try:
        import scanpy as sc

        sc.settings.seed = s
    except Exception:
        pass
    return s


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8", errors="replace")).hexdigest()


def _try_pip_freeze() -> tuple[str | None, str | None]:
    try:
        out = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        if out.returncode == 0 and out.stdout.strip():
            return out.stdout.strip(), _sha256_text(out.stdout.strip())
    except Exception as exc:
        log.debug("pip freeze failed: %s", exc)
    return None, None


def _try_conda_export() -> tuple[str | None, str | None]:
    import os

    if not os.environ.get("CONDA_PREFIX"):
        return None, None
    for cmd in (["conda", "env", "export", "--no-builds"], ["mamba", "env", "export", "--no-builds"]):
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=120, check=False)
            if out.returncode == 0 and out.stdout.strip():
                return out.stdout.strip(), _sha256_text(out.stdout.strip())
        except Exception:
            continue
    return None, None


def compute_environment_fingerprint(*, extra_packages: tuple[str, ...] | None = None) -> dict[str, Any]:
    """Hash full pip/conda env when available; fallback to locked package_versions()."""
    from scagent.export_nb import package_versions

    pip_text, pip_hash = _try_pip_freeze()
    conda_text, conda_hash = _try_conda_export()
    pkg = package_versions(extra=extra_packages or ())
    pkg_json = json.dumps(pkg, sort_keys=True, ensure_ascii=False)
    pkg_hash = _sha256_text(pkg_json)
    sources: list[str] = ["package_versions"]
    parts = [f"package_versions:{pkg_hash}"]
    if pip_hash:
        sources.append("pip_freeze")
        parts.append(f"pip:{pip_hash}")
    if conda_hash:
        sources.append("conda_export")
        parts.append(f"conda:{conda_hash}")
    combined = _sha256_text("|".join(parts))
    return {
        "hash": combined,
        "pip_hash": pip_hash,
        "conda_hash": conda_hash,
        "package_versions_hash": pkg_hash,
        "sources": sources,
        "packages": pkg,
        "python": sys.version.split()[0],
        "pip_freeze_lines": len(pip_text.splitlines()) if pip_text else 0,
        "conda_export_lines": len(conda_text.splitlines()) if conda_text else 0,
    }


def summarize_adata(adata, *, step: str | None = None, path: str | Path | None = None) -> dict[str, Any]:
    """Lightweight AnnData summary for provenance (shape + obs/var columns)."""
    out: dict[str, Any] = {
        "n_obs": int(getattr(adata, "n_obs", 0) or 0),
        "n_vars": int(getattr(adata, "n_vars", 0) or 0),
        "obs_columns": sorted(map(str, getattr(adata, "obs", {}).columns)),
        "var_columns": sorted(map(str, getattr(adata, "var", {}).columns)),
        "layers": sorted(map(str, getattr(adata, "layers", {}).keys())),
        "obsm": sorted(map(str, getattr(adata, "obsm", {}).keys())),
        "backed": bool(getattr(adata, "isbacked", False)),
        "raw": adata.raw is not None,
    }
    if step:
        out["step"] = step
    if path:
        out["path"] = str(path)
    dask_meta = (getattr(adata, "uns", None) or {}).get("scagent_dask")
    if dask_meta:
        out["dask"] = dask_meta
    return out


def summarize_h5ad_path(path: str | Path, *, step: str | None = None) -> dict[str, Any] | None:
    path = Path(path)
    if not path.is_file():
        return None
    try:
        import anndata as ad

        adata = ad.read_h5ad(path, backed="r")
        return summarize_adata(adata, step=step, path=path)
    except Exception as exc:
        log.info("summarize_h5ad_path skip %s: %s", path, exc)
        return {"step": step, "path": str(path), "error": str(exc)}


def build_step_provenance(
    workspace: Path,
    *,
    phase: str,
    data_path: str | None = None,
    metrics: dict | None = None,
) -> list[dict[str, Any]]:
    """Collect input/output AnnData summaries for key pipeline steps."""
    ws = Path(workspace)
    metrics = metrics or {}
    steps: list[dict[str, Any]] = []

    if phase in {"qc", "qc_preprocess.py"} or str(phase).startswith("qc"):
        inp = summarize_h5ad_path(data_path, step="input") if data_path else None
        if data_path and not inp:
            inp = {"step": "input", "path": data_path, "note": "not h5ad or unreadable at manifest time"}
        out = metrics.get("adata_out") or summarize_h5ad_path(ws / "adata_qc.h5ad", step="qc_output")
        steps.append({"step": "qc", "input": inp, "output": out})

    if phase in {"downstream", "cluster_annotate.py"} or "cluster" in str(phase):
        inp = metrics.get("adata_in") or summarize_h5ad_path(ws / "adata_qc.h5ad", step="qc_checkpoint")
        out = metrics.get("adata_out") or summarize_h5ad_path(ws / "adata_processed.h5ad", step="processed_output")
        steps.append({"step": "downstream", "input": inp, "output": out})

    if phase in {"interpret", "interpret_pathways.py"}:
        inp = summarize_h5ad_path(ws / "adata_processed.h5ad", step="interpret_input")
        steps.append({"step": "interpret", "input": inp, "output": None})

    return steps


def merge_manifest_provenance(existing: dict[str, Any], new_steps: list[dict[str, Any]]) -> dict[str, Any]:
    prev = list(existing.get("step_provenance") or [])
    by_step = {s.get("step"): s for s in prev if s.get("step")}
    for s in new_steps:
        if s.get("step"):
            by_step[s["step"]] = s
    existing["step_provenance"] = list(by_step.values())
    return existing


def enrich_run_manifest(path: Path, *, extra: dict[str, Any] | None = None) -> dict[str, Any]:
    """Merge environment + step provenance into workspace/run_manifest.json."""
    from datetime import datetime, timezone

    path = Path(path)
    payload: dict[str, Any] = {}
    if path.is_file():
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            payload = {}
    extra = extra or {}
    seed = int(extra.get("seed") or payload.get("seed") or analysis_params()["seed"])
    payload["seed"] = seed
    payload.setdefault("created_at", datetime.now(timezone.utc).isoformat())
    payload["seed_propagation"] = seed_propagation_record(seed)
    payload["environment"] = compute_environment_fingerprint()
    ws = path.parent
    phase = str(extra.get("phase") or payload.get("phase") or "")
    steps = build_step_provenance(
        ws,
        phase=phase,
        data_path=extra.get("data_path") or payload.get("data_path"),
        metrics=extra.get("metrics") or payload.get("metrics") or {},
    )
    if steps:
        payload = merge_manifest_provenance(payload, steps)
    payload.update({k: v for k, v in extra.items() if k not in {"metrics"}})
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return payload
