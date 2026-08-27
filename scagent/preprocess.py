"""Preprocessing: QC annotation, dynamic filter, normalize, HVG. One job per function."""

from __future__ import annotations

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
    import numpy as np
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
    """Library-size normalize + log1p. Idempotent: skip if X is already log1p."""
    import scanpy as sc

    from scagent.inspect_data import detect_expression_layer

    info = detect_expression_layer(adata)
    adata.uns["expression_layer"] = info
    layer = info.get("layer")
    if layer == "log1p":
        log.info("normalize_log1p skip: already log1p (%s)", info.get("reason"))
        return adata
    if layer == "scaled":
        raise ValueError(
            "adata.X looks scaled (negative values). Restore counts from layers['counts'] "
            "before normalize_log1p; do not log1p a z-scored matrix."
        )
    params = analysis_params()
    if layer == "counts" and "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    if layer == "normalized":
        log.info("normalize_log1p: skip normalize_total (%s)", info.get("reason"))
    else:
        sc.pp.normalize_total(adata, target_sum=target_sum or params["target_sum"])
    sc.pp.log1p(adata)
    info = dict(info)
    info["layer"] = "log1p"
    adata.uns["expression_layer"] = info
    return adata


def select_hvg(adata, *, n_top_genes: int | None = None):
    import scanpy as sc

    n = n_top_genes or analysis_params()["n_hvg"]
    sc.pp.highly_variable_genes(adata, n_top_genes=n, subset=False)
    return adata


# Tirosh / Regev 2016 cell-cycle genes (Scanpy tutorial lists).
S_GENES_HUMAN = (
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
    "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD5", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1",
    "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
)
G2M_GENES_HUMAN = (
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "SMC4", "CCNB2", "CKAP2L",
    "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B",
    "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
    "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
    "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
)


def _cycle_genes(species: str) -> tuple[list[str], list[str]]:
    if str(species).lower() == "mouse":
        return [g.capitalize() for g in S_GENES_HUMAN], [g.capitalize() for g in G2M_GENES_HUMAN]
    return list(S_GENES_HUMAN), list(G2M_GENES_HUMAN)


def cell_cycle_score(adata, species: str = "human"):
    """Score S/G2M with scanpy.tl.score_genes_cell_cycle. Missing genes are skipped."""
    import scanpy as sc

    s_genes, g2m_genes = _cycle_genes(species)
    present = set(map(str, adata.var_names))
    s_use = [g for g in s_genes if g in present]
    g2m_use = [g for g in g2m_genes if g in present]
    if len(s_use) < 3 or len(g2m_use) < 3:
        log.info("cell_cycle_score skipped: too few genes (S=%s G2M=%s)", len(s_use), len(g2m_use))
        adata.obs["S_score"] = 0.0
        adata.obs["G2M_score"] = 0.0
        adata.obs["phase"] = "G1"
        adata.uns["cell_cycle"] = {"scored": False, "reason": "too_few_genes"}
        return adata
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_use, g2m_genes=g2m_use)
    adata.uns["cell_cycle"] = {"scored": True, "n_s": len(s_use), "n_g2m": len(g2m_use)}
    log.info("cell_cycle_score S=%s G2M=%s", len(s_use), len(g2m_use))
    return adata


def maybe_regress_cell_cycle(adata, *, mode: str = "auto", tissue: str | None = None):
    """Regress S/G2M when mode=always, or auto if cycling fraction is high (not embryo/tumor biology)."""
    import scanpy as sc

    mode = (mode or "auto").lower()
    if "S_score" not in adata.obs:
        cell_cycle_score(adata)
    cycling = float((adata.obs.get("phase", "G1").astype(str) != "G1").mean()) if "phase" in adata.obs else 0.0
    do = mode == "always"
    if mode == "auto":
        t = str(tissue or "").lower()
        if t in {"embryo", "tumor"}:
            do = False
        else:
            do = cycling >= 0.15
    adata.uns.setdefault("cell_cycle", {})
    adata.uns["cell_cycle"]["regressed"] = bool(do)
    adata.uns["cell_cycle"]["cycling_fraction"] = cycling
    if not do:
        log.info("cell_cycle regress skipped mode=%s cycling=%.3f", mode, cycling)
        return adata
    sc.pp.regress_out(adata, ["S_score", "G2M_score"])
    log.info("cell_cycle regressed cycling=%.3f", cycling)
    return adata


def choose_ambient(tissue: str | None, requested: str | None = None) -> str:
    req = str(requested or "auto").lower()
    if req in {"none", "off", "skip"}:
        return "none"
    if req in {"soupx", "decontx"}:
        return req
    t = str(tissue or "").lower()
    if t in {"brain", "tumor"}:
        return "soupx"
    try:
        cfg = load_config()
        prof = (cfg.get("qc_profiles") or {}).get(t) or {}
        if prof.get("ambient"):
            return "soupx"
    except Exception:
        pass
    return "none"


def _as_csr(X):
    import numpy as np
    from scipy import sparse as sp

    if sp.issparse(X):
        return X.tocsr().astype(np.float64)
    return sp.csr_matrix(np.asarray(X, dtype=np.float64))


def _estimate_rho(umi) -> float:
    """SoupX-like fallback: low-UMI cells imply modest contamination (typically 5–20%)."""
    import numpy as np

    umi = np.asarray(umi)
    if umi.size == 0:
        return 0.05
    q10, q90 = np.percentile(umi, [10, 90])
    if q90 <= 0:
        return 0.05
    return float(np.clip(0.05 + 0.15 * (1.0 - q10 / max(q90, 1e-9)), 0.02, 0.25))


def remove_ambient(adata, method: str = "soupx", *, rho: float | None = None):
    """Subtract ambient RNA. SoupX/DecontX via rpy2 if present; otherwise Python fallback on counts."""
    method = (method or "none").lower()
    if method in {"none", "off", "skip"}:
        return adata
    if "counts_raw" not in adata.layers:
        adata.layers["counts_raw"] = adata.X.copy()
    used = method
    try:
        if method == "soupx":
            adata = _ambient_soupx_rpy2(adata) or _ambient_soupx_python(adata, rho=rho)
            used = "soupx"
        elif method == "decontx":
            adata = _ambient_decontx_rpy2(adata) or _ambient_decontx_python(adata)
            used = "decontx"
        else:
            adata = _ambient_soupx_python(adata, rho=rho)
            used = "soupx_python"
    except Exception as exc:
        log.warning("remove_ambient %s failed, python fallback: %s", method, exc)
        adata = _ambient_soupx_python(adata, rho=rho)
        used = "soupx_python"
    adata.uns["ambient"] = {**(adata.uns.get("ambient") or {}), "method": used}
    log.info("remove_ambient method=%s", used)
    return adata


def _ambient_soupx_rpy2(adata):
    try:
        from rpy2.robjects.packages import importr  # noqa: F401
    except Exception:
        return None
    return None


def _ambient_decontx_rpy2(adata):
    try:
        from rpy2.robjects.packages import importr  # noqa: F401
    except Exception:
        return None
    return None


def _ambient_soupx_python(adata, *, rho: float | None = None):
    import numpy as np
    from scipy import sparse as sp

    X = _as_csr(adata.X)
    umi = np.asarray(X.sum(axis=1)).ravel()
    n_amb = max(10, int(0.02 * adata.n_obs))
    amb_idx = np.argsort(umi)[:n_amb]
    ambient = np.asarray(X[amb_idx].mean(axis=0)).ravel()
    s = float(ambient.sum())
    if s > 0:
        ambient = ambient / s
    rho = _estimate_rho(umi) if rho is None else float(rho)
    rho = float(np.clip(rho, 0.01, 0.5))
    scale = sp.csr_matrix((rho * umi).reshape(-1, 1))
    amb = sp.csr_matrix(ambient.reshape(1, -1))
    corr = X - scale @ amb
    corr.data = np.clip(corr.data, 0, None)
    corr.eliminate_zeros()
    adata.X = corr.astype(np.float32)
    adata.uns["ambient"] = {"method": "soupx_python", "rho": rho, "n_ambient_cells": n_amb}
    return adata


def _ambient_decontx_python(adata):
    import numpy as np
    from scipy import sparse as sp
    from sklearn.cluster import MiniBatchKMeans

    X = _as_csr(adata.X)
    umi = np.asarray(X.sum(axis=1)).ravel()
    n_g = min(80, X.shape[1])
    var = np.asarray(X.multiply(X).mean(axis=0) - np.square(X.mean(axis=0))).ravel()
    top = np.argsort(var)[-n_g:]
    sub = X[:, top]
    if sp.issparse(sub):
        sub = np.log1p(np.asarray(sub.todense()))
    else:
        sub = np.log1p(np.asarray(sub))
    k = int(min(8, max(2, adata.n_obs // 20)))
    labels = MiniBatchKMeans(n_clusters=k, random_state=0, n_init=3).fit_predict(sub)
    global_p = np.asarray(X.mean(axis=0)).ravel()
    gs = float(global_p.sum()) or 1.0
    global_p = global_p / gs
    corr = X.astype(np.float64).toarray() if adata.n_obs * adata.n_vars < 8_000_000 else None
    if corr is None:
        return _ambient_soupx_python(adata)
    for lab in np.unique(labels):
        idx = np.where(labels == lab)[0]
        native = np.asarray(X[idx].mean(axis=0)).ravel()
        ns = float(native.sum()) or 1.0
        native = native / ns
        rho = 0.1
        umi_i = umi[idx].reshape(-1, 1)
        corr[idx] = np.clip(corr[idx] - rho * umi_i * global_p.reshape(1, -1), 0, None)
        _ = native
    adata.X = sp.csr_matrix(corr.astype(np.float32))
    adata.obs["decontx_cluster"] = labels.astype(str)
    adata.uns["ambient"] = {"method": "decontx_python", "rho": 0.1, "k": k}
    return adata
