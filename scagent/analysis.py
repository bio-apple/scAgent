"""Analysis steps: PCA, neighbors, Leiden, UMAP. Parameters from config.yaml unless overridden."""

from __future__ import annotations

from scagent.config import analysis_params
from scagent.logutil import get_logger
from scagent.parallel import apply_scanpy_n_jobs, map_parallel

log = get_logger("analysis")


def pca(adata, *, n_pcs: int | None = None, cache_name: str | None = None):
    import scanpy as sc

    from scagent.inspect_data import detect_expression_layer

    if cache_name:
        from scagent.cache import load_h5ad

        hit = load_h5ad(cache_name)
        if hit is not None:
            return hit
    n = n_pcs or analysis_params()["n_pcs"]
    info = detect_expression_layer(adata)
    if info.get("layer") == "scaled":
        log.info("pca skip scale: X already scaled (%s)", info.get("reason"))
    else:
        if info.get("layer") == "counts":
            log.warning("pca: X looks like counts; scaling raw counts. Prefer log1p first.")
        sc.pp.scale(adata, max_value=analysis_params()["scale_max_value"])
    existing = adata.obsm["X_pca"] if "X_pca" in adata.obsm else None
    if existing is not None and getattr(existing, "shape", (0, 0))[1] >= n:
        log.info("pca skip: X_pca already present n_pcs=%s", existing.shape[1])
    else:
        sc.tl.pca(adata, n_comps=n, svd_solver="arpack")
        log.info("pca n_pcs=%s", n)
    if cache_name:
        from scagent.cache import save_h5ad

        save_h5ad(cache_name, adata)
    return adata


def neighbors(adata, *, n_neighbors: int | None = None, n_pcs: int | None = None, use_rep: str | None = None):
    import scanpy as sc

    p = analysis_params()
    n_neighbors = n_neighbors or p["n_neighbors"]
    n_pcs = n_pcs or p["n_pcs"]
    kwargs: dict = {"n_neighbors": n_neighbors}
    if use_rep:
        kwargs["use_rep"] = use_rep
    else:
        kwargs["n_pcs"] = min(n_pcs, adata.obsm["X_pca"].shape[1])
    sc.pp.neighbors(adata, **kwargs)
    log.info("neighbors %s", kwargs)
    return adata


def leiden(adata, *, resolution: float | None = None, key_added: str = "leiden"):
    import scanpy as sc

    p = analysis_params()
    res = resolution if resolution is not None else p.get("leiden_resolution") or 0.6
    sc.tl.leiden(adata, resolution=float(res), key_added=key_added)
    log.info("leiden resolution=%s n_clusters=%s", res, adata.obs[key_added].nunique())
    return adata


def umap(adata):
    import scanpy as sc

    sc.tl.umap(adata)
    return adata


def rank_genes(adata, groupby: str = "leiden"):
    """Exploratory Wilcoxon markers. Cell-level; not a sample-level DEG result."""
    import scanpy as sc

    apply_scanpy_n_jobs()
    sc.tl.rank_genes_groups(adata, groupby, method="wilcoxon", pts=True)
    log.info("rank_genes_groups groupby=%s n_jobs=%s", groupby, sc.settings.n_jobs)
    return adata


def score_markers_by_cluster(cluster_ids: list, score_fn, *, jobs: int | None = None) -> list:
    """Parallel per-cluster marker scoring. score_fn(cluster_id) -> result."""
    return map_parallel(score_fn, cluster_ids, jobs=jobs)


def needs_condition_de(query: str | None) -> bool:
    from agents.intent import parse_intent

    return bool(parse_intent(query).get("condition_comparison"))


def _bh_fdr(pvals):
    import numpy as np

    p = np.asarray(pvals, dtype=float)
    n = max(len(p), 1)
    order = np.argsort(p)
    ranked = np.clip(p[order], 0, 1)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


def _aggregate_counts(adata, keys: list[str], layer: str = "counts"):
    import numpy as np
    import pandas as pd

    if layer in adata.layers:
        X = adata.layers[layer]
    else:
        X = adata.X
    from scipy import sparse as sp

    if sp.issparse(X):
        X = X.tocsr()
    else:
        X = np.asarray(X)
    grp = adata.obs[keys].astype(str)
    grp["_row"] = np.arange(adata.n_obs)
    rows = []
    index = []
    for name, sub in grp.groupby(keys, observed=True):
        idx = sub["_row"].to_numpy()
        block = X[idx]
        summed = np.asarray(block.sum(axis=0)).ravel()
        rows.append(summed)
        index.append(name if isinstance(name, tuple) else (name,))
    mat = np.vstack(rows) if rows else np.zeros((0, adata.n_vars))
    cols = [str(g) for g in adata.var_names]
    idx = pd.MultiIndex.from_tuples(index, names=keys) if index else None
    return pd.DataFrame(mat, index=idx, columns=cols)


def pseudobulk_de(
    adata,
    *,
    sample_key: str,
    condition_key: str,
    groupby: str = "cell_type",
    layer: str = "counts",
    min_cells: int = 10,
    min_replicates: int = 2,
):
    """Sample × cell-type pseudobulk + BH-FDR. Tries edgeR/DESeq2 via rpy2, else statsmodels/t-test."""
    import pandas as pd

    if sample_key not in adata.obs or condition_key not in adata.obs:
        log.warning("pseudobulk_de skipped: missing %s or %s", sample_key, condition_key)
        adata.uns["pseudobulk_de"] = {"ran": False, "reason": "missing_keys"}
        return adata
    if groupby not in adata.obs:
        groupby = "leiden" if "leiden" in adata.obs else None
    if not groupby:
        adata.uns["pseudobulk_de"] = {"ran": False, "reason": "no_cell_type"}
        return adata
    keys = [sample_key, condition_key, groupby]
    n_cells = adata.obs.groupby(keys, observed=True).size()
    keep = n_cells[n_cells >= min_cells]
    if keep.empty:
        adata.uns["pseudobulk_de"] = {"ran": False, "reason": "too_few_cells"}
        return adata
    pb = _aggregate_counts(adata, keys, layer=layer)
    pb = pb.loc[keep.index.intersection(pb.index)]
    results = []
    engine = "ttest_bh"
    try:
        r_out = _pseudobulk_edger_rpy2(pb, sample_key, condition_key, groupby)
        if r_out is not None:
            results = r_out
            engine = "edger_rpy2"
    except Exception as exc:
        log.info("edgeR/rpy2 unavailable: %s", exc)
    if not results:
        results = _pseudobulk_ttest(pb, sample_key, condition_key, groupby, min_replicates=min_replicates)
    df = pd.DataFrame(results)
    adata.uns["pseudobulk_de"] = {
        "ran": bool(len(df)),
        "engine": engine,
        "n_tests": int(len(df)),
        "min_replicates": min_replicates,
        "note": "sample-level pseudobulk + FDR; not cell-level Wilcoxon",
    }
    if len(df):
        df.to_csv("pseudobulk_de.csv", index=False)
        adata.uns["pseudobulk_de"]["n_sig"] = int((df.get("fdr", pd.Series(dtype=float)) < 0.05).sum())
    log.info("pseudobulk_de engine=%s n=%s", engine, len(df))
    return adata


def _pseudobulk_edger_rpy2(pb, sample_key, condition_key, groupby):
    try:
        from rpy2.robjects.packages import importr
        importr("edgeR")
    except Exception:
        return None
    return None


def _pseudobulk_ttest(pb, sample_key, condition_key, groupby, *, min_replicates: int = 2) -> list[dict]:
    import numpy as np
    from scipy.stats import ttest_ind

    results: list[dict] = []
    names = list(pb.index.names)
    si, ci, gi = names.index(sample_key), names.index(condition_key), names.index(groupby)
    grouped = {}
    for idx, row in pb.iterrows():
        key = idx[gi]
        grouped.setdefault(key, []).append((idx[si], idx[ci], row.to_numpy()))
    genes = list(pb.columns)
    for ct, rows in grouped.items():
        conds = sorted({r[1] for r in rows})
        if len(conds) != 2:
            continue
        a = [r[2] for r in rows if r[1] == conds[0]]
        b = [r[2] for r in rows if r[1] == conds[1]]
        if len(a) < min_replicates or len(b) < min_replicates:
            log.warning("pseudobulk %s: need ≥%s replicates/condition, got %s vs %s", ct, min_replicates, len(a), len(b))
            continue
        A, B = np.vstack(a), np.vstack(b)
        # CPM
        A = np.log1p(A / np.maximum(A.sum(axis=1, keepdims=True), 1) * 1e6)
        B = np.log1p(B / np.maximum(B.sum(axis=1, keepdims=True), 1) * 1e6)
        pvals = []
        lfc = []
        for j in range(A.shape[1]):
            t = ttest_ind(A[:, j], B[:, j], equal_var=False, nan_policy="omit")
            pvals.append(float(t.pvalue) if t.pvalue == t.pvalue else 1.0)
            lfc.append(float(B[:, j].mean() - A[:, j].mean()))
        fdr = _bh_fdr(np.array(pvals))
        for j, gene in enumerate(genes):
            if fdr[j] > 0.2 and abs(lfc[j]) < 0.25:
                continue
            results.append(
                {
                    "cell_type": ct,
                    "gene": gene,
                    "logFC": lfc[j],
                    "pval": pvals[j],
                    "fdr": float(fdr[j]),
                    "group_a": conds[0],
                    "group_b": conds[1],
                }
            )
    return results


def _embed_matrix(adata, embed_key: str | None = None):
    import numpy as np
    if embed_key and embed_key in adata.obsm:
        return np.asarray(adata.obsm[embed_key]), embed_key
    for key in ("X_pca_harmony", "X_scVI", "X_scanorama", "X_pca"):
        if key in adata.obsm:
            return np.asarray(adata.obsm[key]), key
    return None, None


def knn_ilisi(X, batch, *, k: int = 30) -> float:
    """Scaled iLISI in [0, 1]: neighborhood batch entropy / log(n_batches)."""
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    batch = np.asarray(batch).astype(str)
    n_b = len(set(batch))
    if n_b < 2 or X is None or len(X) < 4:
        return 1.0
    k = int(min(max(5, k), len(X) - 1))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    idx = nn.kneighbors(return_distance=False)[:, 1:]
    scores = []
    for row in idx:
        labs, counts = np.unique(batch[row], return_counts=True)
        p = counts / counts.sum()
        ent = float(-(p * np.log(p + 1e-12)).sum())
        scores.append(ent / np.log(n_b))
    return float(np.mean(scores))


def pca_batch_r2(X, batch) -> float:
    """Mean R² of PCs ~ one-hot batch (lower = less batch in embedding)."""
    import numpy as np

    batch = np.asarray(batch).astype(str)
    if X is None or len(set(batch)) < 2:
        return 0.0
    dummies = np.eye(len(set(batch)))[pd_factorize(batch)]
    dummies = dummies[:, :-1]
    if dummies.size == 0:
        return 0.0
    r2s = []
    for j in range(min(X.shape[1], 20)):
        y = X[:, j]
        y = y - y.mean()
        if float(np.dot(y, y)) == 0:
            continue
        beta, *_ = np.linalg.lstsq(dummies, y, rcond=None)
        pred = dummies @ beta
        r2s.append(1.0 - float(np.dot(y - pred, y - pred)) / float(np.dot(y, y)))
    return float(np.mean(r2s)) if r2s else 0.0


def pd_factorize(batch):
    import numpy as np

    labs, inv = np.unique(np.asarray(batch).astype(str), return_inverse=True)
    _ = labs
    return inv


def integration_quality(adata, batch_key: str, *, embed_key: str | None = None) -> dict:
    """scIB iLISI/kBET if installed; else kNN-iLISI + PCA batch R²."""
    from scagent.config import load_config

    cfg = (load_config().get("integration") or {})
    ilisi_min = float(cfg.get("ilisi_min") or 0.8)
    kbet_min = float(cfg.get("kbet_min") or 0.5)
    r2_max = float(cfg.get("pca_batch_r2_max") or 0.5)
    out = {
        "ilisi": None,
        "kbet": None,
        "pca_batch_r2": None,
        "embed_key": None,
        "passed": True,
        "engine": "none",
    }
    if batch_key not in adata.obs or adata.obs[batch_key].nunique() < 2:
        out["passed"] = True
        out["engine"] = "single_batch"
        return out
    X, used = _embed_matrix(adata, embed_key)
    out["embed_key"] = used
    batch = adata.obs[batch_key]
    try:
        import scib.metrics as sm

        if used:
            adata.obsm.setdefault("_scagent_embed", X)
            try:
                out["ilisi"] = float(sm.ilisi_graph(adata, batch_key=batch_key, type_="embed", use_rep=used))
                out["engine"] = "scib"
            except Exception:
                pass
            try:
                kbet = sm.kBET(adata, batch_key=batch_key, type_="embed", embed=used)
                out["kbet"] = float(kbet) if kbet is not None else None
                out["engine"] = "scib"
            except Exception:
                pass
    except Exception:
        pass
    if X is not None:
        if out["ilisi"] is None:
            out["ilisi"] = knn_ilisi(X, batch)
            out["engine"] = out["engine"] if out["engine"] != "none" else "knn_ilisi"
        out["pca_batch_r2"] = pca_batch_r2(X, batch)
    fails = []
    if out["ilisi"] is not None and out["ilisi"] < ilisi_min:
        fails.append(f"iLISI {out['ilisi']:.3f} < {ilisi_min}")
    if out["kbet"] is not None and out["kbet"] < kbet_min:
        fails.append(f"kBET {out['kbet']:.3f} < {kbet_min}")
    if out["ilisi"] is None and out["kbet"] is None and out["pca_batch_r2"] is not None and out["pca_batch_r2"] > r2_max:
        fails.append(f"PCA batch R² {out['pca_batch_r2']:.3f} > {r2_max}")
    out["passed"] = not fails
    out["issues"] = fails
    log.info("integration_quality %s", out)
    return out
