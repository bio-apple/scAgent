"""Analysis steps: PCA, neighbors, Leiden, UMAP. Parameters from config.yaml unless overridden."""

from __future__ import annotations

from scagent.config import analysis_params
from scagent.logutil import get_logger
from scagent.parallel import apply_scanpy_n_jobs, map_parallel

log = get_logger("analysis")


def pca(adata, *, n_pcs: int | None = None, cache_name: str | None = None):
    import scanpy as sc

    if cache_name:
        from scagent.cache import load_h5ad

        hit = load_h5ad(cache_name)
        if hit is not None:
            return hit
    n = n_pcs or analysis_params()["n_pcs"]
    sc.pp.scale(adata, max_value=analysis_params()["scale_max_value"])
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
