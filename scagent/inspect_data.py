from __future__ import annotations

import json
from pathlib import Path
from typing import Any

# log1p(CP10K/CPM) almost never exceeds this; raw UMI max usually does.
_LOG1P_MAX = 20.0
_INT_TOL = 1e-6


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


_CC_GENES = {
    "MKI67",
    "PCNA",
    "TOP2A",
    "CDK1",
    "CCNB1",
    "CCNB2",
    "CCNA2",
    "MCM2",
    "MCM6",
    "TYMS",
    "HMGB2",
    "UBE2C",
}


def gene_composition(genes: list[str]) -> dict[str, Any]:
    """Presence of MT / ribo / HB / cell-cycle gene symbols (not expression)."""
    sample = [str(g) for g in genes[: min(len(genes), 8000)]]
    n = max(len(sample), 1)
    up = [g.upper() for g in sample]
    n_mt = sum(g.startswith("MT-") or g.startswith("MT.") for g in up)
    n_ribo = sum(g.startswith(("RPS", "RPL")) for g in up)
    n_hb = sum(g.startswith("HB") and not g.startswith("HBP") for g in up)
    n_cc = sum(g in _CC_GENES for g in up)
    return {
        "n_mt_genes": n_mt,
        "n_ribo_genes": n_ribo,
        "n_hb_genes": n_hb,
        "n_cc_genes": n_cc,
        "frac_ribo_symbols": round(n_ribo / n, 4),
        "need_hb_qc": n_hb >= 3,
        "need_cell_cycle": n_cc >= 2,
    }


def _sample_x_stats(X, *, max_cells: int = 256) -> dict[str, Any]:
    """Max/min/sparsity/integer-likeness on a cell subset (works with backed/sparse)."""
    import numpy as np

    if X is None or getattr(X, "shape", (0, 0))[0] == 0:
        return {"x_min": None, "x_max": None, "sparsity": None, "is_integer_like": None, "n_sampled": 0}
    n_obs = int(X.shape[0])
    take = min(n_obs, max_cells)
    sub = X[:take]
    try:
        from scipy import sparse as sp

        is_sp = sp.issparse(sub)
    except Exception:
        is_sp = False
    if is_sp:
        sub = sub.tocsr()
        size = int(sub.shape[0]) * int(sub.shape[1])
        nnz = int(sub.nnz)
        data = np.asarray(sub.data, dtype=np.float64) if nnz else np.zeros(1, dtype=np.float64)
        x_max = float(data.max()) if nnz else 0.0
        x_min = float(data.min()) if nnz else 0.0
        if nnz < size:
            x_min = min(x_min, 0.0)
        sparsity = 1.0 - (nnz / max(size, 1))
        frac_int = float(np.mean(np.abs(data - np.round(data)) < _INT_TOL)) if nnz else 1.0
    else:
        arr = np.asarray(sub, dtype=np.float64)
        x_min = float(np.min(arr))
        x_max = float(np.max(arr))
        sparsity = float(np.mean(arr == 0))
        frac_int = float(np.mean(np.abs(arr - np.round(arr)) < _INT_TOL))
    return {
        "x_min": x_min,
        "x_max": x_max,
        "sparsity": sparsity,
        "is_integer_like": bool(frac_int >= 0.9),
        "n_sampled": take,
    }


def _uns_normalization_flags(uns: Any) -> tuple[bool, bool]:
    if uns is None:
        return False, False
    try:
        keys = {str(k).lower() for k in uns.keys()}
    except Exception:
        return False, False
    has_log1p = "log1p" in keys
    has_norm = bool(keys & {"normalize_total", "normalization", "n_counts_norm"})
    return has_log1p, has_norm


def detect_expression_layer(adata, *, max_cells: int = 256) -> dict[str, Any]:
    """Classify adata.X as counts / normalized / log1p / scaled.

    Priority: uns records (scanpy log1p / normalize_total), then negatives (scaled),
    then max + integer-likeness + sparsity.
    """
    uns = getattr(adata, "uns", None)
    has_log1p, has_norm = _uns_normalization_flags(uns)
    layers: list[str] = []
    try:
        layers = [str(k) for k in (getattr(adata, "layers", None) or {}).keys()]
    except Exception:
        layers = []
    stats = _sample_x_stats(getattr(adata, "X", None), max_cells=max_cells)
    x_min, x_max = stats["x_min"], stats["x_max"]
    sparsity = stats["sparsity"]
    integer_like = stats["is_integer_like"]
    has_neg = x_min is not None and x_min < -0.01
    dense = sparsity is not None and sparsity < 0.3

    if has_neg or (dense and x_min is not None and x_max is not None and x_min < 0):
        layer, reason = "scaled", "X has negative values (z-score / scale)"
    elif has_log1p:
        layer, reason = "log1p", "adata.uns contains log1p"
    elif has_norm and not has_log1p:
        layer, reason = "normalized", "adata.uns records size-factor normalization without log1p"
    elif x_max is not None and x_max > _LOG1P_MAX:
        layer, reason = "counts", f"X max={x_max:.3g} > {_LOG1P_MAX} (raw/CPM-like)"
    elif integer_like:
        layer, reason = "counts", "X is integer-like (UMI/counts)"
    elif x_max is not None and x_max <= _LOG1P_MAX:
        layer, reason = "log1p", f"X max={x_max:.3g} ≤ {_LOG1P_MAX} and not integer-like"
    else:
        layer, reason = "unknown", "could not classify X"

    return {
        "layer": layer,
        "reason": reason,
        "x_min": x_min,
        "x_max": x_max,
        "sparsity": sparsity,
        "is_integer_like": integer_like,
        "uns_log1p": has_log1p,
        "uns_normalized": has_norm,
        "layers": layers,
        "n_sampled": stats["n_sampled"],
    }


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
            genes = [str(g) for g in adata.var_names[:8000]]
            if species is None:
                meta["species"] = _detect_species_from_genes(genes)
            comp = gene_composition(genes)
            meta.update(comp)
            if comp.get("need_hb_qc"):
                meta["notes"].append("检测到血红蛋白基因，QC 建议纳入 hb。")
            if comp.get("need_cell_cycle"):
                meta["notes"].append("检测到细胞周期基因，建议 score_genes_cell_cycle；regress 由 config.qc.regress_cell_cycle 决定。")
            if int(meta.get("n_cells") or 0) >= 100_000:
                meta["notes"].append("细胞数 ≥10 万：模板将走 backed/CSR/n_jobs。")
            preferred = None
            if extra:
                preferred = extra.get("sample_key") or extra.get("batch_key")
            keys = [k for k in (preferred, "sample", "batch", "donor", "orig.ident", "orig_ident") if k]
            for key in keys:
                if key in adata.obs.columns:
                    n = int(adata.obs[key].nunique())
                    meta["sample_key"] = key
                    meta["n_samples"] = n
                    meta["need_batch_correction"] = n > 1
                    break
            for ckey in ("condition", "group", "treatment", "status", "disease"):
                if ckey in adata.obs.columns:
                    meta["condition_key"] = ckey
                    break
            sk, ck = meta.get("sample_key"), meta.get("condition_key")
            if sk and ck and sk in adata.obs.columns and ck in adata.obs.columns:
                try:
                    pairs = adata.obs[[sk, ck]].drop_duplicates()
                    n_s = int(adata.obs[sk].nunique())
                    n_c = int(adata.obs[ck].nunique())
                    if n_s > 1 and n_s == n_c and len(pairs) == n_s:
                        meta["batch_condition_confounded"] = True
                        meta["notes"].append(
                            "样本与条件 1:1 共线：auto 将跳过整合，以免把处理效应当批次抹掉。"
                        )
                except Exception:
                    pass
            try:
                expr = detect_expression_layer(adata)
                meta["expression_layer"] = expr.get("layer")
                meta["x_max"] = expr.get("x_max")
                meta["x_min"] = expr.get("x_min")
                meta["sparsity"] = expr.get("sparsity")
                meta["uns_log1p"] = expr.get("uns_log1p")
                meta["uns_normalized"] = expr.get("uns_normalized")
                meta["layers"] = expr.get("layers") or []
                layer = expr.get("layer")
                if layer == "log1p":
                    meta["notes"].append(
                        "X 已是 log1p，将跳过 normalize_total/log1p，避免重复归一化。"
                    )
                elif layer == "scaled":
                    meta["notes"].append(
                        "X 已是 scaled（含负值）。不要再 scale/log1p；请从 layers['counts'] 恢复。"
                    )
                elif layer == "normalized":
                    meta["notes"].append("X 已做 size-factor 归一化但未 log1p；将只补 log1p。")
            except Exception as exc:
                meta["notes"].append(f"表达层检测失败: {exc}")
            try:
                from scagent.config import performance_params

                meta["backed_recommended"] = int(meta.get("n_cells") or 0) >= int(
                    performance_params()["backed_threshold_cells"]
                )
            except Exception:
                meta["backed_recommended"] = False
            try:
                adata.file.close()
            except Exception:
                pass
        except ImportError:
            meta["notes"].append("未安装 anndata，跳过 h5ad 内部检查。")
        except Exception as exc:  # pragma: no cover - IO edge
            meta["notes"].append(f"读取 h5ad 失败: {exc}")

    elif path.suffix.lower() in {".rds", ".h5seurat"}:
        meta["notes"].append(
            "检测到 Seurat 对象。Python 路径用 scagent.io.read_single_cell 转为 AnnData（需 R + zellkonverter 或 rpy2）。"
        )

    if meta["n_samples"] > 1:
        meta["need_batch_correction"] = True
    return meta


detect_species_from_genes = _detect_species_from_genes
detect_platform = _detect_platform


def metadata_json(meta: dict[str, Any]) -> str:
    return json.dumps(meta, ensure_ascii=False, indent=2)
