from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from scagent.io import (
    count_tsv_rows,
    discover_samples,
    parse_data_spec,
    peek_10x_h5_shape,
    peek_loom,
    resolve_10x_h5,
    resolve_10x_matrix_dir,
    sample_label,
)

# log1p(CP10K/CPM) almost never exceeds this; raw UMI max usually does.
_LOG1P_MAX = 20.0
_INT_TOL = 1e-6

BATCH_KEYS = (
    "sample",
    "batch",
    "donor",
    "orig.ident",
    "orig_ident",
    "library_id",
    "sample_id",
    "Batch",
    "SAMPLE",
)


def _detect_platform(path: Path) -> str:
    name = path.name.lower()
    if "parse" in name:
        return "parse"
    if name.endswith(".loom"):
        return "loom"
    if name.endswith(".rds") or name.endswith(".h5seurat"):
        return "seurat"
    if resolve_10x_matrix_dir(path) is not None or resolve_10x_h5(path) is not None:
        return "10x"
    if path.is_dir() and (
        (path / "matrix.mtx").exists()
        or (path / "matrix.mtx.gz").exists()
        or (path / "filtered_feature_bc_matrix").exists()
        or (path / "outs" / "filtered_feature_bc_matrix").exists()
        or (path / "outs" / "filtered_feature_bc_matrix.h5").exists()
    ):
        return "10x"
    if name.endswith(".h5") or "10x" in name or name.endswith(".h5ad"):
        return "10x"
    if path.is_dir() and len(discover_samples(path)) > 1:
        return "multi"
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


def _preferred_batch_keys(extra: dict[str, Any] | None) -> list[str]:
    preferred = None
    if extra:
        preferred = extra.get("sample_key") or extra.get("batch_key")
    keys = [k for k in (preferred, *BATCH_KEYS) if k]
    seen: set[str] = set()
    out: list[str] = []
    for k in keys:
        if k in seen:
            continue
        seen.add(k)
        out.append(k)
    return out


def _apply_batch_from_columns(
    meta: dict[str, Any],
    columns: list[str],
    nunique: dict[str, int],
    extra: dict[str, Any] | None,
) -> None:
    for key in _preferred_batch_keys(extra):
        if key not in columns:
            continue
        n = int(nunique.get(key) or 0)
        if n <= 0:
            continue
        meta["sample_key"] = key
        meta["n_samples"] = n
        meta["need_batch_correction"] = n > 1
        break
    for ckey in ("condition", "group", "treatment", "status", "disease"):
        if ckey in columns:
            meta["condition_key"] = ckey
            break


def _append_batch_auto_note(meta: dict[str, Any]) -> None:
    if int(meta.get("n_samples") or 1) > 1:
        meta["need_batch_correction"] = True
    if meta.get("batch_condition_confounded"):
        return
    if not meta.get("need_batch_correction"):
        return
    if any("将触发批次校正" in str(n) for n in meta.get("notes") or []):
        return
    key = meta.get("sample_key") or "sample"
    n = meta.get("n_samples")
    meta["notes"].append(
        f"检测到批次列 {key!r}（n_samples={n}），auto 将触发批次校正"
        "（小数据 Harmony；≥10 万细胞或 ≥8 样本 scVI；也可指定 Scanorama/BBKNN）。"
    )


def _inspect_h5ad(path: Path, meta: dict[str, Any], *, species: str | None, extra: dict[str, Any] | None) -> None:
    try:
        import anndata as ad  # type: ignore
    except ImportError:
        meta["notes"].append("未安装 anndata，跳过 h5ad 内部检查。")
        return
    try:
        adata = ad.read_h5ad(path, backed="r")
        meta["n_cells"] = int(adata.n_obs)
        meta["n_genes"] = int(adata.n_vars)
        meta["obs_columns"] = list(map(str, adata.obs.columns))
        genes = [str(g) for g in adata.var_names[:8000]]
        if species is None or meta.get("species") in {None, "unknown"}:
            meta["species"] = _detect_species_from_genes(genes)
        comp = gene_composition(genes)
        meta.update(comp)
        if comp.get("need_hb_qc"):
            meta["notes"].append("检测到血红蛋白基因，QC 建议纳入 hb。")
        if comp.get("need_cell_cycle"):
            meta["notes"].append(
                "检测到细胞周期基因，建议 score_genes_cell_cycle；regress 由 config.qc.regress_cell_cycle 决定。"
            )
        if int(meta.get("n_cells") or 0) >= 100_000:
            meta["notes"].append("细胞数 ≥10 万：模板将走 backed/CSR/n_jobs。")
        nunique = {c: int(adata.obs[c].nunique()) for c in adata.obs.columns}
        _apply_batch_from_columns(meta, meta["obs_columns"], nunique, extra)
        sk, ck = meta.get("sample_key"), meta.get("condition_key")
        if sk and ck and sk in adata.obs.columns and ck in adata.obs.columns:
            try:
                pairs = adata.obs[[sk, ck]].drop_duplicates()
                n_s = int(adata.obs[sk].nunique())
                n_c = int(adata.obs[ck].nunique())
                if n_s > 1 and n_s == n_c and len(pairs) == n_s:
                    meta["batch_condition_confounded"] = True
                    meta["notes"].append("样本与条件 1:1 共线：auto 将跳过整合，以免把处理效应当批次抹掉。")
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
                meta["notes"].append("X 已是 log1p，将跳过 normalize_total/log1p，避免重复归一化。")
            elif layer == "scaled":
                meta["notes"].append("X 已是 scaled（含负值）。不要再 scale/log1p；请从 layers['counts'] 恢复。")
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
    except Exception as exc:  # pragma: no cover - IO edge
        meta["notes"].append(f"读取 h5ad 失败: {exc}")


def _inspect_one_path(path: Path, meta: dict[str, Any], *, species: str | None, extra: dict[str, Any] | None) -> None:
    if meta.get("platform") in {None, "unknown"}:
        meta["platform"] = _detect_platform(path)
    suf = path.suffix.lower()
    if suf == ".h5ad":
        _inspect_h5ad(path, meta, species=species, extra=extra)
        return
    if suf == ".loom":
        peek = peek_loom(path)
        if peek.get("n_cells"):
            meta["n_cells"] = peek["n_cells"]
        if peek.get("n_genes"):
            meta["n_genes"] = peek["n_genes"]
        cols = list(peek.get("obs_columns") or [])
        meta["obs_columns"] = cols
        _apply_batch_from_columns(meta, cols, dict(peek.get("obs_nunique") or {}), extra)
        genes = list(peek.get("genes") or [])
        if genes:
            if species is None or meta.get("species") in {None, "unknown"}:
                meta["species"] = _detect_species_from_genes(genes)
            meta.update(gene_composition(genes))
        meta["notes"].append("检测到 Loom。执行时用 anndata.read_loom（需 loompy）。")
        return
    if suf in {".rds", ".h5seurat"}:
        meta["platform"] = "seurat"
        meta["notes"].append(
            "检测到 Seurat 对象。Python 路径用 scagent.io.read_single_cell 转为 AnnData（需 R + zellkonverter 或 rpy2）。"
        )
        return
    mtx = resolve_10x_matrix_dir(path)
    h5 = resolve_10x_h5(path)
    if mtx is not None or h5 is not None:
        meta["platform"] = "10x"
        meta["expression_layer"] = "counts"
        if mtx is not None:
            n_bc = count_tsv_rows(mtx, ("barcodes.tsv.gz", "barcodes.tsv"))
            n_ft = count_tsv_rows(mtx, ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))
            if n_bc:
                meta["n_cells"] = n_bc
            if n_ft:
                meta["n_genes"] = n_ft
            meta["notes"].append(f"Cell Ranger / 10x mtx: {mtx}")
        elif h5 is not None:
            n_cells, n_genes = peek_10x_h5_shape(h5)
            if n_cells:
                meta["n_cells"] = n_cells
            if n_genes:
                meta["n_genes"] = n_genes
            meta["notes"].append(f"Cell Ranger 10x h5: {h5}")
        meta["notes"].append("10x 视为 raw counts；执行时 sc.read_10x_mtx / read_10x_h5。")
        return
    if suf == ".h5":
        meta["platform"] = "10x"
        n_cells, n_genes = peek_10x_h5_shape(path)
        if n_cells:
            meta["n_cells"] = n_cells
        if n_genes:
            meta["n_genes"] = n_genes


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
        meta["notes"] = list(meta.get("notes") or [])

    if not data_path:
        meta["notes"].append("未提供数据路径，仅根据任务描述规划。")
        return meta

    specs = parse_data_spec(data_path)
    missing = [p for p in specs if not p.exists()]
    if not specs or (missing and len(specs) == len(missing)):
        path = Path(data_path)
        meta["exists"] = False
        meta["notes"].append(f"路径不存在: {missing[0] if missing else path}")
        if platform:
            meta["platform"] = platform
        return meta
    if missing:
        meta["notes"].append("部分路径不存在: " + ", ".join(str(p) for p in missing))
    specs = [p for p in specs if p.exists()]

    samples: list[Path] = []
    seen: set[str] = set()
    for spec in specs:
        for p in discover_samples(spec) or [spec]:
            key = str(p.resolve()) if p.exists() else str(p)
            if key in seen:
                continue
            seen.add(key)
            samples.append(p)
    if not samples:
        meta["exists"] = False
        meta["notes"].append(f"路径不存在: {data_path}")
        return meta

    meta["exists"] = True
    meta["input_paths"] = [str(p) for p in samples]
    multi = len(samples) > 1
    if multi:
        user_key = None
        if extra:
            user_key = extra.get("sample_key") or extra.get("batch_key")
        meta["sample_key"] = user_key or "sample"
        meta["n_samples"] = len(samples)
        meta["need_batch_correction"] = True
        labels = [sample_label(p) for p in samples]
        shown = ", ".join(labels[:8]) + ("…" if len(labels) > 8 else "")
        meta["notes"].append(f"将拼接 {len(samples)} 个输入（{shown}）到 obs[{meta['sample_key']!r}]。")
        meta["platform"] = meta["platform"] if meta.get("platform") not in {None, "unknown"} else "multi"

    first = samples[0]
    if meta.get("platform") in {None, "unknown"}:
        meta["platform"] = _detect_platform(first)
    _inspect_one_path(first, meta, species=species, extra=extra)

    if multi:
        user_key = None
        if extra:
            user_key = extra.get("sample_key") or extra.get("batch_key")
        meta["sample_key"] = user_key or "sample"
        meta["n_samples"] = len(samples)
        meta["need_batch_correction"] = True
        total = 0
        any_cells = False
        for sample in samples:
            piece: dict[str, Any] = {"notes": [], "platform": "unknown"}
            _inspect_one_path(sample, piece, species=species, extra=extra)
            if piece.get("n_cells"):
                total += int(piece["n_cells"])
                any_cells = True
        if any_cells:
            meta["n_cells"] = total
        # first path already inspected; extra notes from siblings except duplicates
        # (first was inspected twice — acceptable for inspect)

    try:
        from scagent.trajectory import inspect_trajectory_hints

        hints = inspect_trajectory_hints(
            tissue=meta.get("tissue"),
            layers=list(meta.get("layers") or []),
            n_obs=meta.get("n_cells"),
            query=(extra or {}).get("task") or (extra or {}).get("query"),
        )
        meta.update({k: v for k, v in hints.items() if k != "trajectory_notes"})
        meta["notes"].extend(hints.get("trajectory_notes") or [])
    except Exception:
        pass

    _append_batch_auto_note(meta)
    return meta


detect_species_from_genes = _detect_species_from_genes
detect_platform = _detect_platform


def metadata_json(meta: dict[str, Any]) -> str:
    return json.dumps(meta, ensure_ascii=False, indent=2)
