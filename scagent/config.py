from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml
from dotenv import load_dotenv

from scagent.logutil import setup_logging

REPO_ROOT = Path(__file__).resolve().parents[1]
_CACHE: dict[str, Any] | None = None


def load_config(path: Path | None = None, *, reload: bool = False) -> dict[str, Any]:
    """Load YAML config. Path: argument, else SCAGENT_CONFIG, else repo config.yaml.
    API keys stay in environment variables, never in the YAML file."""
    global _CACHE
    if _CACHE is not None and not reload and path is None:
        return _CACHE
    load_dotenv(REPO_ROOT / ".env")
    import os

    cfg_path = path or os.getenv("SCAGENT_CONFIG") or (REPO_ROOT / "config.yaml")
    cfg_path = Path(cfg_path)
    with cfg_path.open(encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    cfg["_root"] = str(REPO_ROOT)
    cfg["_config_path"] = str(cfg_path)
    log_cfg = cfg.get("logging") or {}
    setup_logging(level=log_cfg.get("level"), log_file=_abs(cfg, log_cfg.get("file")))
    _CACHE = cfg
    return cfg


def cfg_get(cfg: dict[str, Any] | None, dotted: str, default: Any = None) -> Any:
    data: Any = cfg if cfg is not None else load_config()
    for part in dotted.split("."):
        if not isinstance(data, dict) or part not in data:
            return default
        data = data[part]
    return default if data is None else data


def analysis_params(cfg: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = cfg or load_config()
    p = cfg.get("params") or {}
    a = cfg.get("analysis") or {}
    return {
        "n_pcs": int(p.get("n_pcs") or 40),
        "n_neighbors": int(p.get("n_neighbors") or 15),
        "n_hvg": int(p.get("n_hvg") or 2000),
        "hvg_flavor": str(p.get("hvg_flavor") or "seurat_v3"),
        "leiden_resolution": p.get("leiden_resolution"),
        "leiden_resolutions": list(p.get("leiden_resolutions") or [0.2, 0.4, 0.6, 0.8, 1.0]),
        "target_sum": float(p.get("target_sum") or 1e4),
        "scale_max_value": float(p.get("scale_max_value") or 10),
        "seed": int(a.get("seed") if a.get("seed") is not None else p.get("seed") or 0),
    }


def performance_params(cfg: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = cfg or load_config()
    p = cfg.get("performance") or {}
    n_jobs = p.get("n_jobs")
    return {
        "n_jobs": int(n_jobs) if n_jobs is not None else -1,
        "backed_threshold_cells": int(p.get("backed_threshold_cells") or 250000),
        "cache": True if p.get("cache") is None else bool(p.get("cache")),
    }


def resolve_path(cfg: dict[str, Any], key: str) -> Path:
    rel = cfg["paths"][key]
    p = Path(rel)
    if p.is_absolute():
        return p
    return REPO_ROOT / p


def _abs(cfg: dict[str, Any], rel: str | None) -> Path | None:
    if not rel:
        return None
    p = Path(rel)
    if p.is_absolute():
        return p
    return Path(cfg.get("_root") or REPO_ROOT) / p
