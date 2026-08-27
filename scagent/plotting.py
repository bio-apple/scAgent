"""Visualization helpers. Figure directory from config paths.workspace/figures unless overridden."""

from __future__ import annotations

from pathlib import Path

from scagent.config import load_config, resolve_path
from scagent.logutil import get_logger

log = get_logger("plotting")


def _figdir(figdir: str | Path | None) -> Path:
    if figdir:
        p = Path(figdir)
    else:
        cfg = load_config()
        p = resolve_path(cfg, "workspace") / "figures"
    p.mkdir(parents=True, exist_ok=True)
    return p


def apply_figdir(figdir: str | Path | None = None):
    import scanpy as sc

    d = _figdir(figdir)
    sc.settings.figdir = str(d)
    return d


def qc_violin(adata, *, save: str = "_qc_violin.png", figdir: str | Path | None = None):
    import scanpy as sc

    apply_figdir(figdir)
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=save,
        show=False,
    )
    log.info("qc violin %s", save)


def qc_scatter(adata, *, figdir: str | Path | None = None):
    import scanpy as sc

    apply_figdir(figdir)
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save="_qc_scatter_counts.png", show=False)
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save="_qc_scatter_mt.png", show=False)
    log.info("qc scatter")


def umap(adata, color, *, save: str = "_umap.png", figdir: str | Path | None = None):
    import scanpy as sc

    apply_figdir(figdir)
    sc.pl.umap(adata, color=color, save=save, show=False)
    log.info("umap %s", save)
