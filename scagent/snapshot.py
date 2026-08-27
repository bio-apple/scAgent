"""Lazy AnnData snapshots: hardlink full matrices, store obs/obsm deltas when X is unchanged.

AnnData is never kept in AgentState. Branching forks the on-disk index, not RAM copies of X.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scagent.config import REPO_ROOT, load_config
from scagent.logutil import get_logger

log = get_logger("snapshot")

_INDEX = "index.json"


def _cache_dir(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    rel = (cfg.get("paths") or {}).get("cache") or ".cache"
    p = Path(rel)
    if not p.is_absolute():
        p = Path(cfg.get("_root") or REPO_ROOT) / p
    p.mkdir(parents=True, exist_ok=True)
    return p


def snapshots_root(cfg: dict | None = None, thread_id: str | None = None) -> Path:
    root = _cache_dir(cfg) / "snapshots"
    if thread_id:
        root = root / str(thread_id)
    root.mkdir(parents=True, exist_ok=True)
    return root


def hardlink_or_copy(src: Path, dest: Path) -> str:
    """Share inode when possible so a branch does not duplicate a multi-GB h5ad."""
    src = Path(src)
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not src.is_file():
        raise FileNotFoundError(src)
    if dest.exists():
        try:
            if dest.samefile(src):
                return "link"
        except OSError:
            pass
        dest.unlink()
    try:
        os.link(src, dest)
        return "link"
    except OSError:
        shutil.copy2(src, dest)
        return "copy"


def x_fingerprint(path: Path) -> str | None:
    """Identify X without loading the matrix into RAM (h5py slice of the data vector)."""
    path = Path(path)
    try:
        import h5py
        import numpy as np
    except ImportError:
        return f"size:{path.stat().st_size}"
    try:
        with h5py.File(path, "r") as f:
            x = f.get("X")
            if x is None:
                return None
            if hasattr(x, "keys") and "data" in x:
                data = x["data"]
                shape = tuple(int(s) for s in (x.attrs.get("shape") or data.shape))
                n = min(512, int(data.shape[0]))
                head = np.asarray(data[:n]).tobytes() if n else b""
            else:
                shape = tuple(int(s) for s in x.shape)
                n0 = min(8, int(shape[0]) if shape else 0)
                n1 = min(32, int(shape[1]) if len(shape) > 1 else 0)
                head = np.asarray(x[:n0, :n1]).tobytes() if n0 and n1 else b""
            return f"{shape}:{hashlib.sha256(head).hexdigest()[:16]}"
    except Exception as exc:
        log.debug("x_fingerprint failed: %s", exc)
        return f"size:{path.stat().st_size}"


def _load_index(thread_id: str, cfg: dict | None = None) -> dict[str, Any]:
    path = snapshots_root(cfg, thread_id) / _INDEX
    if not path.is_file():
        return {"thread_id": thread_id, "entries": []}
    return json.loads(path.read_text(encoding="utf-8"))


def _save_index(payload: dict[str, Any], cfg: dict | None = None) -> Path:
    tid = payload["thread_id"]
    path = snapshots_root(cfg, tid) / _INDEX
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


def list_snapshots(thread_id: str, cfg: dict | None = None) -> list[dict[str, Any]]:
    return list(_load_index(thread_id, cfg).get("entries") or [])


def _write_obs_delta(src: Path, dest_dir: Path) -> dict[str, str]:
    """Persist obs (+ obsm keys) without touching X. Requires anndata; falls back to empty."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    files: dict[str, str] = {}
    try:
        import anndata as ad
    except ImportError:
        return files
    adata = ad.read_h5ad(src, backed="r")
    obs_path = dest_dir / "obs.csv.gz"
    adata.obs.to_csv(obs_path, compression="gzip")
    files["obs"] = str(obs_path)
    try:
        import numpy as np

        obsm = {k: np.asarray(v) for k, v in dict(adata.obsm).items()}
        if obsm:
            npz = dest_dir / "obsm.npz"
            np.savez_compressed(npz, **obsm)
            files["obsm"] = str(npz)
    except Exception as exc:
        log.debug("obsm delta skipped: %s", exc)
    return files


def record_h5ad(
    src: Path,
    *,
    step: str,
    thread_id: str,
    parent: dict[str, Any] | None = None,
    params: dict[str, Any] | None = None,
    cfg: dict | None = None,
) -> dict[str, Any]:
    src = Path(src)
    root = snapshots_root(cfg, thread_id)
    step_dir = root / step
    step_dir.mkdir(parents=True, exist_ok=True)
    fp = x_fingerprint(src)
    parent_fp = (parent or {}).get("x_fingerprint")
    parent_path = (parent or {}).get("path")
    kind = "link"
    dest_h5ad = step_dir / src.name
    extra: dict[str, str] = {}
    if parent_path and parent_fp and fp and fp == parent_fp and Path(parent_path).is_file():
        extra = _write_obs_delta(src, step_dir / "delta")
        if extra.get("obs"):
            kind = "delta"
            dest_h5ad = Path(parent_path)
    if kind != "delta":
        kind = hardlink_or_copy(src, dest_h5ad)
    entry = {
        "id": f"{thread_id}:{step}",
        "thread_id": thread_id,
        "step": step,
        "kind": kind,
        "path": str(dest_h5ad),
        "source": str(src),
        "parent": (parent or {}).get("id"),
        "parent_path": str(parent_path) if parent_path else None,
        "x_fingerprint": fp,
        "bytes": int(src.stat().st_size),
        "stored_bytes": int(Path(extra["obs"]).stat().st_size) if kind == "delta" else int(dest_h5ad.stat().st_size),
        "params": params or {},
        "created_at": datetime.now(timezone.utc).isoformat(),
        **extra,
    }
    idx = _load_index(thread_id, cfg)
    idx["entries"] = [e for e in idx.get("entries") or [] if e.get("step") != step]
    idx["entries"].append(entry)
    _save_index(idx, cfg)
    log.info("snapshot %s kind=%s step=%s stored=%s", thread_id, kind, step, entry["stored_bytes"])
    return entry


def record_phase(
    workspace: Path,
    phase: str,
    *,
    thread_id: str | None,
    params: dict[str, Any] | None = None,
    cfg: dict | None = None,
) -> list[dict[str, Any]]:
    """Replace full-file copies of phase h5ads with hardlinks / obs deltas."""
    workspace = Path(workspace)
    tid = str(thread_id or "default")
    names = {"qc": ("adata_qc.h5ad",), "downstream": ("adata_processed.h5ad",)}
    wanted = names.get(phase) or ("adata_qc.h5ad", "adata_processed.h5ad")
    existing = list_snapshots(tid, cfg)
    parent = None
    for e in reversed(existing):
        if e.get("step") == "qc" or e.get("kind") in {"link", "copy", "full"}:
            parent = e
            if e.get("step") == "qc":
                break
    out: list[dict[str, Any]] = []
    for name in wanted:
        src = workspace / name
        if not src.is_file():
            continue
        step = "qc" if "qc" in name else phase
        use_parent = parent if step != "qc" else None
        out.append(record_h5ad(src, step=step, thread_id=tid, parent=use_parent, params=params, cfg=cfg))
        if step == "qc":
            parent = out[-1]
    return out


def open_lazy(entry: dict[str, Any], *, backed: bool = True):
    """Open matrix with backed='r'. Delta snapshots keep X on the parent file; obs overlay is separate."""
    from scagent.io import read_h5ad

    path = entry.get("parent_path") if entry.get("kind") == "delta" else entry.get("path")
    adata = read_h5ad(path, backed=backed)
    overlay = None
    obs_path = entry.get("obs")
    if obs_path and Path(obs_path).is_file():
        import pandas as pd

        overlay = pd.read_csv(obs_path, index_col=0, compression="gzip")
    return {"adata": adata, "obs_overlay": overlay, "kind": entry.get("kind")}


def fork_branch(
    src_thread: str,
    dest_thread: str,
    *,
    from_step: str | None = None,
    cfg: dict | None = None,
) -> dict[str, Any]:
    """New thread_id sharing parent h5ad inodes. Does not copy X into RAM or duplicate GB files."""
    src_entries = list_snapshots(src_thread, cfg)
    if from_step:
        keep: list[dict[str, Any]] = []
        for e in src_entries:
            keep.append(e)
            if e.get("step") == from_step:
                break
        src_entries = keep
    dest_root = snapshots_root(cfg, dest_thread)
    new_entries: list[dict[str, Any]] = []
    for e in src_entries:
        step = e["step"]
        cloned = dict(e)
        cloned["id"] = f"{dest_thread}:{step}"
        cloned["thread_id"] = dest_thread
        cloned["forked_from"] = e.get("id")
        if e.get("kind") == "delta":
            delta_src = Path(e["obs"]).parent if e.get("obs") else None
            if delta_src and delta_src.is_dir():
                dest_delta = dest_root / step / "delta"
                if dest_delta.exists():
                    shutil.rmtree(dest_delta)
                shutil.copytree(delta_src, dest_delta)
                cloned["obs"] = str(dest_delta / Path(e["obs"]).name) if e.get("obs") else None
                if e.get("obsm"):
                    cloned["obsm"] = str(dest_delta / Path(e["obsm"]).name)
            cloned["path"] = e.get("path")
            cloned["parent_path"] = e.get("parent_path") or e.get("path")
        else:
            src_file = Path(e["path"])
            dest_file = dest_root / step / src_file.name
            if src_file.is_file():
                cloned["kind"] = hardlink_or_copy(src_file, dest_file)
                cloned["path"] = str(dest_file)
        new_entries.append(cloned)
    payload = {"thread_id": dest_thread, "parent_thread": src_thread, "from_step": from_step, "entries": new_entries}
    _save_index(payload, cfg)
    log.info("forked branch %s → %s steps=%s", src_thread, dest_thread, [e["step"] for e in new_entries])
    return payload


def checkout(thread_id: str, step: str, workspace: Path, cfg: dict | None = None) -> Path | None:
    """Point workspace at a snapshot via hardlink (no extra GB copy)."""
    for e in list_snapshots(thread_id, cfg):
        if e.get("step") != step:
            continue
        src = Path(e["parent_path"] if e.get("kind") == "delta" else e["path"])
        if not src.is_file():
            return None
        name = "adata_qc.h5ad" if step == "qc" else "adata_processed.h5ad"
        dest = Path(workspace) / name
        hardlink_or_copy(src, dest)
        return dest
    return None
