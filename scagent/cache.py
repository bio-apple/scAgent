"""Disk cache for AnnData checkpoints, JSON, and LLM responses."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from scagent.config import load_config, performance_params, resolve_path
from scagent.logutil import get_logger

log = get_logger("cache")


def cache_enabled(cfg: dict | None = None) -> bool:
    return bool(performance_params(cfg)["cache"])


def cache_dir(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    try:
        d = resolve_path(cfg, "cache")
    except (KeyError, TypeError):
        d = Path(cfg.get("_root") or ".") / ".cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


def key_hash(*parts: Any) -> str:
    blob = "||".join(str(p) for p in parts).encode("utf-8", errors="replace")
    return hashlib.sha256(blob).hexdigest()[:24]


def save_json(name: str, payload: Any, cfg: dict | None = None) -> Path | None:
    if not cache_enabled(cfg):
        return None
    path = cache_dir(cfg) / f"{name}.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    log.debug("cache write %s", path)
    return path


def load_json(name: str, cfg: dict | None = None) -> Any | None:
    if not cache_enabled(cfg):
        return None
    path = cache_dir(cfg) / f"{name}.json"
    if not path.exists():
        return None
    log.info("cache hit %s", path)
    return json.loads(path.read_text(encoding="utf-8"))


def save_h5ad(name: str, adata, cfg: dict | None = None) -> Path | None:
    if not cache_enabled(cfg):
        return None
    path = cache_dir(cfg) / f"{name}.h5ad"
    adata.write_h5ad(path)
    log.info("cache write %s", path)
    return path


def load_h5ad(name: str, cfg: dict | None = None):
    if not cache_enabled(cfg):
        return None
    path = cache_dir(cfg) / f"{name}.h5ad"
    if not path.exists():
        return None
    import anndata as ad

    log.info("cache hit %s", path)
    return ad.read_h5ad(path)


def llm_key(messages) -> str:
    parts = []
    for m in messages or []:
        role = getattr(m, "type", None) or getattr(m, "role", "") or type(m).__name__
        content = getattr(m, "content", m)
        parts.append(f"{role}:{content}")
    return "llm_" + key_hash(*parts)
