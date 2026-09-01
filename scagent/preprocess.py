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


def filter_dynamic(
    adata,
    *,
    method: str | None = None,
    nmads: int = 5,
    percentile: dict | None = None,
    sample_key: str | None = None,
    per_sample: bool | None = None,
):
    """MAD and/or percentile filter. No default mito%<5.

    When ``per_sample`` is True (default if sample_key has >1 levels), thresholds
    are computed within each sample then OR-merged — avoids one sample dominating MAD.
    """
    import numpy as np
    from scipy.stats import median_abs_deviation

    cfg = load_config()
    method = method or (cfg.get("qc") or {}).get("method") or "mad"
    percentile = percentile or (cfg.get("qc") or {}).get("percentile") or {}
    n_before = adata.n_obs
    sk = sample_key
    if per_sample is None:
        per_sample = bool(sk and sk in adata.obs.columns and int(adata.obs[sk].nunique()) > 1)
    if per_sample and sk and sk in adata.obs.columns:
        outlier = np.zeros(n_before, dtype=bool)
        for sample in adata.obs[sk].astype(str).unique():
            cell_mask = (adata.obs[sk].astype(str) == str(sample)).to_numpy()
            sub = adata[cell_mask]
            outlier[cell_mask] = _outlier_mask(sub, method=method, nmads=nmads, percentile=percentile)
        adata.obs["outlier"] = outlier
        out = adata[~outlier].copy()
        log.info(
            "filter_dynamic method=%s per_sample=%s removed=%s remaining=%s",
            method,
            sk,
            n_before - out.n_obs,
            out.n_obs,
        )
        return out
    outlier = _outlier_mask(adata, method=method, nmads=nmads, percentile=percentile)
    adata.obs["outlier"] = outlier
    out = adata[~outlier].copy()
    log.info("filter_dynamic method=%s removed=%s remaining=%s", method, n_before - out.n_obs, out.n_obs)
    return out


def _outlier_mask(adata, *, method: str, nmads: int, percentile: dict):
    import numpy as np
    from scipy.stats import median_abs_deviation

    n = adata.n_obs
    x_mt = adata.obs["pct_counts_mt"].to_numpy()
    outlier = np.zeros(n, dtype=bool)

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
    return outlier


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


def normalize_expression(adata, *, method: str = "log1p", target_sum: float | None = None):
    """Normalize counts. method: log1p | pearson | sctransform.

    ``sctransform`` tries Seurat SCTransform via R; on failure falls back to Pearson
    residuals and records honesty metadata (never silently claims SCT success).
    """
    method = str(method or "log1p").lower().strip()
    if method in {"none", "off", "skip"}:
        adata.uns["normalization"] = {"method": "none", "applied": False}
        return adata
    if method in {"log1p", "lognorm", "lognormalize"}:
        normalize_log1p(adata, target_sum=target_sum)
        adata.uns["normalization"] = {"method": "log1p", "applied": True}
        return adata
    if method in {"pearson", "pearson_residuals"}:
        return _normalize_pearson(adata)
    if method in {"sctransform", "sct"}:
        ok = _normalize_sctransform_r(adata)
        if ok:
            return adata
        log.warning("SCTransform unavailable; falling back to pearson_residuals")
        print("SCAGENT_WARN: SCTransform unavailable; used pearson_residuals (not Seurat SCT)")
        out = _normalize_pearson(adata)
        info = dict(out.uns.get("normalization") or {})
        info.update({"requested": "sctransform", "fallback": "pearson_residuals"})
        out.uns["normalization"] = info
        return out
    raise ValueError(f"unknown normalization method: {method}")


def _normalize_pearson(adata):
    import scanpy as sc

    from scagent.inspect_data import detect_expression_layer

    info = detect_expression_layer(adata)
    if "counts" not in adata.layers and info.get("layer") == "counts":
        adata.layers["counts"] = adata.X.copy()
    layer = "counts" if "counts" in adata.layers else None
    try:
        sc.experimental.pp.normalize_pearson_residuals(adata, layer=layer)
    except Exception as exc:
        log.warning("pearson residuals failed (%s); falling back to log1p", exc)
        normalize_log1p(adata)
        adata.uns["normalization"] = {
            "method": "log1p",
            "applied": True,
            "requested": "pearson",
            "note": str(exc),
        }
        return adata
    adata.uns["normalization"] = {"method": "pearson_residuals", "applied": True}
    adata.uns["expression_layer"] = {"layer": "scaled", "reason": "pearson_residuals"}
    return adata


def _normalize_sctransform_r(adata) -> bool:
    """Attempt Seurat::SCTransform via Rscript. Returns True on success."""
    import shutil
    import subprocess
    import tempfile
    from pathlib import Path

    if not shutil.which("Rscript"):
        return False
    r_script = Path(__file__).resolve().parent / "r" / "sctransform.R"
    if not r_script.is_file():
        return False
    tmp = Path(tempfile.mkdtemp(prefix="scagent_sct_"))
    try:
        inp = tmp / "in.h5ad"
        out = tmp / "out.h5ad"
        adata.write_h5ad(inp)
        proc = subprocess.run(
            ["Rscript", str(r_script), str(inp), str(out)],
            capture_output=True,
            text=True,
            timeout=600,
        )
        if proc.returncode != 0 or not out.exists():
            log.debug("sctransform.R failed: %s", proc.stderr or proc.stdout)
            return False
        import anndata as ad

        loaded = ad.read_h5ad(out)
        adata.X = loaded.X
        for col in loaded.obs.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = loaded.obs[col].to_numpy()
        if "counts" not in adata.layers and "counts" in (loaded.layers or {}):
            adata.layers["counts"] = loaded.layers["counts"]
        adata.uns["normalization"] = {"method": "sctransform", "applied": True, "backend": "seurat_r"}
        adata.uns["expression_layer"] = {"layer": "scaled", "reason": "sctransform"}
        return True
    except Exception as exc:
        log.debug("sctransform failed: %s", exc)
        return False
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def select_hvg(adata, *, n_top_genes: int | None = None, batch_key: str | None = None, flavor: str | None = None, random_state: int | None = None):
    """HVG on raw counts with seurat_v3 when possible (Heumos 2023). Multi-batch uses batch_key union."""
    import scanpy as sc

    p = analysis_params()
    n = n_top_genes or p["n_hvg"]
    rs = int(random_state if random_state is not None else p["seed"])
    flavor = flavor or p.get("hvg_flavor") or "seurat_v3"
    kwargs: dict = {"n_top_genes": n, "subset": False, "flavor": flavor}
    if flavor == "seurat_v3" and "counts" in adata.layers:
        kwargs["layer"] = "counts"
    elif flavor == "seurat_v3" and "counts" not in adata.layers:
        from scagent.inspect_data import detect_expression_layer

        if detect_expression_layer(adata).get("layer") != "counts":
            kwargs["flavor"] = "seurat"
            log.info("select_hvg: no counts layer; fallback flavor=seurat on X")
    if batch_key and batch_key in adata.obs and adata.obs[batch_key].nunique() > 1:
        kwargs["batch_key"] = batch_key
    try:
        sc.pp.highly_variable_genes(adata, **kwargs)
    except Exception as exc:
        log.warning("select_hvg %s failed (%s); fallback seurat", kwargs.get("flavor"), exc)
        fb = {"n_top_genes": n, "subset": False, "flavor": "seurat"}
        if batch_key and batch_key in adata.obs and adata.obs[batch_key].nunique() > 1:
            fb["batch_key"] = batch_key
        sc.pp.highly_variable_genes(adata, **fb)
        kwargs = fb
    adata.uns["scagent_hvg"] = {"flavor": kwargs.get("flavor"), "n_top_genes": n, "random_state": rs}
    log.info("select_hvg flavor=%s n=%s batch_key=%s random_state=%s", kwargs.get("flavor"), n, kwargs.get("batch_key"), rs)
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


def ambient_backend_available(method: str = "soupx") -> bool:
    """True only when a real SoupX/DecontX R binding can run (currently unwired → False)."""
    primary = "decontx" if "decontx" in str(method).lower() else "soupx"
    try:
        from rpy2.robjects.packages import importr

        if primary == "soupx":
            importr("SoupX")
            # Import alone is insufficient until adjustCounts is wired.
            return False
        importr("celda")
        return False
    except Exception:
        return False


def choose_ambient(tissue: str | None, requested: str | None = None) -> str:
    """Pick ambient method. Auto never requests SoupX/DecontX until a real backend is wired."""
    req = str(requested or "auto").lower()
    if req in {"none", "off", "skip"}:
        return "none"
    if req in {"soupx_heuristic", "decontx_heuristic"}:
        return req
    if req in {"soupx", "decontx"}:
        # Explicit user request is preserved; remove_ambient reports applied=False if unavailable.
        return req
    # auto: do not advertise SoupX for brain/tumor until R binding exists
    if ambient_backend_available("soupx"):
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


def remove_ambient(
    adata,
    method: str = "soupx",
    *,
    rho: float | None = None,
    allow_heuristic: bool = False,
):
    """Subtract ambient RNA via real SoupX/DecontX when available.

    Heuristic Python fallback does **not** mutate counts unless ``allow_heuristic=True``
    or method is explicitly ``soupx_heuristic`` / ``decontx_heuristic``.
    """
    method = (method or "none").lower()
    if method in {"none", "off", "skip"}:
        adata.uns["ambient"] = {"method": "none", "applied": False}
        return adata
    want_heuristic = method in {"soupx_heuristic", "decontx_heuristic"} or allow_heuristic
    primary = "decontx" if "decontx" in method else "soupx"
    if "counts_raw" not in adata.layers:
        adata.layers["counts_raw"] = adata.X.copy()
    try:
        real = _ambient_soupx_rpy2(adata) if primary == "soupx" else _ambient_decontx_rpy2(adata)
        if real is not None:
            adata = real
            info = dict(adata.uns.get("ambient") or {})
            info.update({"method": primary, "applied": True, "backend": "rpy2"})
            adata.uns["ambient"] = info
            log.info("remove_ambient method=%s applied=True", primary)
            return adata
    except Exception as exc:
        log.warning("remove_ambient %s rpy2 failed: %s", primary, exc)
    if want_heuristic:
        if primary == "decontx":
            adata = _ambient_decontx_python(adata)
        else:
            adata = _ambient_soupx_python(adata, rho=rho)
        info = dict(adata.uns.get("ambient") or {})
        info["applied"] = True
        info["requested"] = method
        adata.uns["ambient"] = info
        log.warning("remove_ambient applied heuristic fallback method=%s", info.get("method"))
        return adata
    adata.uns["ambient"] = {
        "method": f"{primary}_unavailable",
        "requested": method,
        "applied": False,
        "note": "SoupX/DecontX not available; counts unchanged (set allow_heuristic=True to force heuristic)",
    }
    log.warning("remove_ambient skipped: %s unavailable; counts unchanged", primary)
    return adata


def _ambient_soupx_rpy2(adata):
    """Real SoupX via rpy2. Returns None until wired to SoupX::autoEstCont / adjustCounts."""
    try:
        from rpy2.robjects.packages import importr  # noqa: F401
    except Exception:
        return None
    # Package import alone is not SoupX — do not claim success.
    return None


def _ambient_decontx_rpy2(adata):
    """Real DecontX via rpy2. Returns None until wired to celda::decontX."""
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
    adata.uns["ambient"] = {
        "method": "soupx_heuristic",
        "rho": rho,
        "n_ambient_cells": n_amb,
        "applied": True,
    }
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
    adata.uns["ambient"] = {"method": "decontx_heuristic", "rho": 0.1, "k": k, "applied": True}
    return adata
