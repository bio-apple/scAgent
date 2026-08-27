"""Preprocessing: QC annotation, dynamic filter, normalize, HVG. One job per function."""

from __future__ import annotations

import numpy as np

from scagent.config import analysis_params, load_config
from scagent.logutil import get_logger

log = get_logger("preprocess")


def annotate_qc_genes(adata, *, species: str = "human"):
    mt = "MT-" if species != "mouse" else "mt-"
    adata.var["mt"] = adata.var_names.str.startswith(mt)
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]", case=False, regex=True)
    return adata


def calculate_qc(adata, *, qc_vars: list[str] | None = None):
    import scanpy as sc

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=qc_vars or ["mt"], percent_top=None, log1p=True, inplace=True
    )
    return adata


def filter_dynamic(adata, *, method: str | None = None, nmads: int = 5, percentile: dict | None = None):
    """MAD and/or percentile filter. No default mito%<5."""
    from scipy.stats import median_abs_deviation

    cfg = load_config()
    method = method or (cfg.get("qc") or {}).get("method") or "mad"
    percentile = percentile or (cfg.get("qc") or {}).get("percentile") or {}
    n_before = adata.n_obs
    x_mt = adata.obs["pct_counts_mt"].to_numpy()
    outlier = np.zeros(n_before, dtype=bool)

    def mad_mask(values, side="two"):
        med = np.median(values)
        mad = median_abs_deviation(values)
        if mad == 0:
            return np.zeros(len(values), dtype=bool)
        if side == "high":
            return values > med + nmads * mad
        return (values < med - nmads * mad) | (values > med + nmads * mad)

    if method in {"mad", "hybrid"}:
        outlier |= mad_mask(x_mt, side="high")
        outlier |= mad_mask(adata.obs["log1p_total_counts"].to_numpy(), side="two")
        outlier |= mad_mask(adata.obs["log1p_n_genes_by_counts"].to_numpy(), side="two")
    if method in {"percentile", "hybrid"}:
        g = adata.obs["n_genes_by_counts"].to_numpy()
        c = adata.obs["total_counts"].to_numpy()
        outlier |= g < np.percentile(g, percentile.get("n_genes_low") or 2)
        outlier |= g > np.percentile(g, percentile.get("n_genes_high") or 98)
        outlier |= c < np.percentile(c, percentile.get("total_counts_low") or 2)
        outlier |= c > np.percentile(c, percentile.get("total_counts_high") or 98)
        outlier |= x_mt > np.percentile(x_mt, percentile.get("pct_mt_high") or 98)
    adata = adata[~outlier].copy()
    log.info("filter_dynamic method=%s removed=%s remaining=%s", method, n_before - adata.n_obs, adata.n_obs)
    return adata


def normalize_log1p(adata, *, target_sum: float | None = None):
    import scanpy as sc

    params = analysis_params()
    sc.pp.normalize_total(adata, target_sum=target_sum or params["target_sum"])
    sc.pp.log1p(adata)
    return adata


def select_hvg(adata, *, n_top_genes: int | None = None):
    import scanpy as sc

    n = n_top_genes or analysis_params()["n_hvg"]
    sc.pp.highly_variable_genes(adata, n_top_genes=n, subset=False)
    return adata
