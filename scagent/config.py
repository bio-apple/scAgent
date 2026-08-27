from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml
from dotenv import load_dotenv

REPO_ROOT = Path(__file__).resolve().parents[1]


def load_config(path: Path | None = None) -> dict[str, Any]:
    load_dotenv(REPO_ROOT / ".env")
    cfg_path = path or REPO_ROOT / "config.yaml"
    with cfg_path.open(encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    cfg["_root"] = str(REPO_ROOT)
    return cfg


def resolve_path(cfg: dict[str, Any], key: str) -> Path:
    rel = cfg["paths"][key]
    p = Path(rel)
    if p.is_absolute():
        return p
    return REPO_ROOT / p
