from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _detect_platform(path: Path) -> str:
    name = path.name.lower()
    if "parse" in name:
        return "parse"
    if path.is_dir() and (
        (path / "matrix.mtx").exists()
        or (path / "matrix.mtx.gz").exists()
        or (path / "filtered_feature_bc_matrix").exists()
    ):
        return "10x"
    if name.endswith(".h5") or "10x" in name or name.endswith(".h5ad"):
        return "10x"
    if name.endswith(".rds") or name.endswith(".h5seurat"):
        return "seurat"
    return "unknown"


def _detect_species_from_genes(genes: list[str]) -> str:
    sample = genes[: min(len(genes), 5000)]
    mt_human = sum(g.startswith("MT-") for g in sample)
    mt_mouse = sum(g.startswith("mt-") for g in sample)
    if mt_human > mt_mouse and mt_human > 0:
        return "human"
    if mt_mouse > mt_human and mt_mouse > 0:
        return "mouse"
    return "unknown"


def inspect_data(
    data_path: str | None,
    *,
    tissue: str | None = None,
    species: str | None = None,
    platform: str | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Read lightweight metadata. Missing scanpy/anndata is not fatal."""
    meta: dict[str, Any] = {
        "data_path": data_path,
        "exists": False,
        "species": species or "unknown",
        "platform": platform or "unknown",
        "tissue": (tissue or "default").lower(),
        "n_cells": None,
        "n_genes": None,
        "n_samples": 1,
        "sample_key": None,
        "obs_columns": [],
        "need_batch_correction": False,
        "notes": [],
    }
    if extra:
        meta.update(extra)

    if not data_path:
        meta["notes"].append("未提供数据路径，仅根据任务描述规划。")
        return meta

    path = Path(data_path)
    meta["exists"] = path.exists()
    if not path.exists():
        meta["notes"].append(f"路径不存在: {path}")
        if platform:
            meta["platform"] = platform
        return meta

    if meta["platform"] == "unknown":
        meta["platform"] = _detect_platform(path)

    if path.suffix.lower() == ".h5ad":
        try:
            import anndata as ad  # type: ignore

            adata = ad.read_h5ad(path, backed="r")
            meta["n_cells"] = int(adata.n_obs)
            meta["n_genes"] = int(adata.n_vars)
            meta["obs_columns"] = list(map(str, adata.obs.columns))
            genes = [str(g) for g in adata.var_names[:5000]]
            if species is None:
                meta["species"] = _detect_species_from_genes(genes)
            for key in ("sample", "batch", "donor", "orig.ident", "orig_ident"):
                if key in adata.obs.columns:
                    n = int(adata.obs[key].nunique())
                    meta["sample_key"] = key
                    meta["n_samples"] = n
                    meta["need_batch_correction"] = n > 1
                    break
            try:
                adata.file.close()
            except Exception:
                pass
        except ImportError:
            meta["notes"].append("未安装 anndata，跳过 h5ad 内部检查。")
        except Exception as exc:  # pragma: no cover - IO edge
            meta["notes"].append(f"读取 h5ad 失败: {exc}")

    if meta["n_samples"] > 1:
        meta["need_batch_correction"] = True
    return meta


def metadata_json(meta: dict[str, Any]) -> str:
    return json.dumps(meta, ensure_ascii=False, indent=2)
