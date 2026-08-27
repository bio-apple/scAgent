"""Large-scale paths: backed I/O, experimental Dask, GPU (scVI / RAPIDS)."""

from __future__ import annotations

from typing import Any

from scagent.config import load_config
from scagent.logutil import get_logger

log = get_logger("performance")


def dask_params(cfg: dict | None = None) -> dict[str, Any]:
    cfg = cfg or load_config()
    p = (cfg.get("performance") or {}).get("dask") or {}
    return {
        "enabled": bool(p.get("enabled")),
        "threshold_cells": int(p.get("threshold_cells") or 500_000),
        "chunk_size": int(p.get("chunk_size") or 10_000),
    }


def gpu_params(cfg: dict | None = None) -> dict[str, Any]:
    cfg = cfg or load_config()
    p = (cfg.get("performance") or {}).get("gpu") or {}
    return {
        "enabled": bool(p.get("enabled")),
        "scvi": True if p.get("scvi") is None else bool(p.get("scvi")),
        "rapids": bool(p.get("rapids")),
    }


def cuda_available() -> bool:
    try:
        import torch

        return bool(torch.cuda.is_available())
    except Exception:
        return False


def scvi_train_kwargs(cfg: dict | None = None) -> dict[str, Any]:
    """Extra kwargs for scvi.model.SCVI.train when GPU is enabled and available."""
    gp = gpu_params(cfg)
    if not gp["enabled"] or not gp["scvi"]:
        return {}
    if not cuda_available():
        log.info("scVI GPU requested but CUDA unavailable; CPU train")
        return {}
    # scvi-tools ≥1.0 uses Lightning accelerator=gpu
    try:
        import scvi

        parts = [p for p in str(getattr(scvi, "__version__", "0")).split(".")[:2] if p.isdigit()]
        ver = tuple(int(x) for x in parts) if parts else (0, 0)
        if ver >= (1, 0):
            return {"accelerator": "gpu", "devices": 1}
    except Exception:
        pass
    return {"use_gpu": True}


def scvi_train_suffix(cfg: dict | None = None) -> str:
    kw = scvi_train_kwargs(cfg)
    if not kw:
        return ""
    return ", " + ", ".join(f"{k}={v!r}" for k, v in kw.items())


def configure_scanpy_dask(adata, cfg: dict | None = None) -> None:
    """Tag AnnData for Scanpy experimental out-of-core / Dask path."""
    dp = dask_params(cfg)
    if not dp["enabled"]:
        return
    n = int(getattr(adata, "n_obs", 0) or 0)
    if n < dp["threshold_cells"]:
        return
    adata.uns["scagent_dask"] = {
        "enabled": True,
        "chunk_size": dp["chunk_size"],
        "experimental": True,
        "note": "Scanpy out-of-core; backed/dask chunking — subset before scale/PCA when needed",
    }
    try:
        import scanpy as sc

        if hasattr(sc.settings, "dask"):
            sc.settings.dask = True
    except Exception:
        pass
    log.info("dask experimental path enabled n_obs=%s chunk=%s", n, dp["chunk_size"])


def try_rapids_neighbors(adata, *, n_neighbors: int, n_pcs: int, use_rep: str | None = None) -> bool:
    """Run RAPIDS/cuML neighbors when configured. Returns True if used."""
    gp = gpu_params()
    if not gp["enabled"] or not gp["rapids"]:
        return False
    if not cuda_available():
        log.info("RAPIDS requested but CUDA unavailable")
        return False
    try:
        import rapids_singlecell as rsc
    except ImportError:
        try:
            import rapids_singlecell.pp as rsc_pp  # type: ignore

            rsc = rsc_pp
        except ImportError:
            log.info("rapids-singlecell not installed")
            return False
    rep = use_rep or "X_pca"
    if rep not in adata.obsm:
        return False
    try:
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=rep)
        adata.uns["scagent_gpu"] = {"rapids_neighbors": True, "use_rep": rep}
        log.info("RAPIDS neighbors n=%s rep=%s", n_neighbors, rep)
        return True
    except Exception as exc:
        log.warning("RAPIDS neighbors failed: %s", exc)
        return False


def try_rapids_umap(adata) -> bool:
    gp = gpu_params()
    if not gp["enabled"] or not gp["rapids"]:
        return False
    if not cuda_available():
        return False
    try:
        import rapids_singlecell as rsc
    except ImportError:
        return False
    try:
        rsc.pp.umap(adata)
        adata.uns.setdefault("scagent_gpu", {})["rapids_umap"] = True
        log.info("RAPIDS umap")
        return True
    except Exception as exc:
        log.warning("RAPIDS umap failed: %s", exc)
        return False
