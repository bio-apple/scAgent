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


def _subsample_idx(n: int, max_cells: int, seed: int = 0):
    import numpy as np

    if n <= max_cells:
        return np.arange(n)
    return np.random.default_rng(seed).choice(n, size=max_cells, replace=False)


def _scatter_batch(xy, labels, path: Path, title: str) -> Path | None:
    import numpy as np

    xy = np.asarray(xy)
    if xy.ndim != 2 or xy.shape[1] < 2 or xy.shape[0] == 0:
        return None
    labels = np.asarray(labels).astype(str)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.2, 4.4), dpi=120)
    labs = list(dict.fromkeys(labels.tolist()))
    cmap = plt.get_cmap("tab20")
    for i, lab in enumerate(labs):
        m = labels == lab
        ax.scatter(xy[m, 0], xy[m, 1], s=6, alpha=0.75, color=cmap(i % 20), label=lab, linewidths=0)
    ax.set_title(title, fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
    if len(labs) <= 12:
        ax.legend(markerscale=2, fontsize=8, frameon=False, loc="best")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    log.info("integration diagnostic %s", path.name)
    return path


def marker_heatmap(
    adata,
    *,
    groupby: str = "leiden",
    n_genes: int = 8,
    save: str = "_marker_heatmap.png",
    figdir: str | Path | None = None,
) -> str | None:
    """Top cluster marker heatmap (exploratory). Requires rank_genes_groups in adata.uns."""
    import scanpy as sc

    if "rank_genes_groups" not in (adata.uns or {}):
        return None
    if groupby not in adata.obs:
        return None
    d = apply_figdir(figdir)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=n_genes,
        groupby=groupby,
        show_gene_labels=True,
        show=False,
    )
    out = d / "marker_heatmap.png"
    plt.savefig(out, bbox_inches="tight", dpi=120)
    plt.close("all")
    log.info("marker heatmap %s", out.name)
    return str(out) if out.is_file() else None


def volcano_from_de_csv(
    csv_path: str | Path,
    *,
    save: str = "volcano.png",
    figdir: str | Path | None = None,
    fdr_cutoff: float = 0.05,
) -> str | None:
    """Volcano plot from pseudobulk_de.csv (gene, logFC, pval, fdr)."""
    import csv

    path = Path(csv_path)
    if not path.is_file():
        return None
    rows: list[dict] = []
    with path.open(encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        for rec in reader:
            try:
                lfc = float(rec.get("logFC") or rec.get("logfc") or "nan")
                pval = float(rec.get("pval") or rec.get("p_value") or "nan")
                fdr = float(rec.get("fdr") or rec.get("padj") or pval)
            except (TypeError, ValueError):
                continue
            if lfc != lfc or pval != pval:
                continue
            rows.append({"gene": rec.get("gene") or "", "logFC": lfc, "pval": max(pval, 1e-300), "fdr": fdr})
    if not rows:
        return None
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    out_dir = _figdir(figdir)
    out = out_dir / save
    x = np.array([r["logFC"] for r in rows])
    y = -np.log10(np.array([r["pval"] for r in rows]))
    sig = np.array([r["fdr"] < fdr_cutoff for r in rows])
    fig, ax = plt.subplots(figsize=(5.5, 4.8), dpi=120)
    ax.scatter(x[~sig], y[~sig], s=8, alpha=0.45, c="#94a3b8", linewidths=0)
    ax.scatter(x[sig], y[sig], s=10, alpha=0.75, c="#dc2626", linewidths=0)
    ax.axhline(-np.log10(0.05), color="#64748b", lw=0.8, ls="--")
    ax.set_xlabel("log2 FC")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Pseudobulk DEG volcano")
    top = sorted(rows, key=lambda r: r["fdr"])[:8]
    for r in top:
        if r["fdr"] >= fdr_cutoff:
            break
        idx = rows.index(r)
        ax.annotate(str(r["gene"]), (x[idx], y[idx]), fontsize=7, alpha=0.85)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    log.info("volcano plot %s", out.name)
    return str(out)


def pathway_bubble_plot(
    terms: list[dict],
    *,
    save: str = "pathway_bubble.png",
    figdir: str | Path | None = None,
    top_n: int = 15,
    fdr_cutoff: float = 0.05,
) -> str | None:
    """Bubble plot: y=pathway, x=-log10(p), size=overlap."""
    if not terms:
        return None
    ranked = sorted(terms, key=lambda t: float(t.get("fdr") or t.get("pval") or 1.0))[:top_n]
    if not ranked:
        return None
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    out_dir = _figdir(figdir)
    out = out_dir / save
    labels = [str(t.get("term") or t.get("name") or "?")[:48] for t in ranked]
    pvals = np.array([max(float(t.get("pval") or 1.0), 1e-300) for t in ranked])
    sizes = np.array([max(int(t.get("overlap") or 1), 1) for t in ranked]) * 18.0
    sig = np.array([float(t.get("fdr") or t.get("pval") or 1.0) < fdr_cutoff for t in ranked])
    y = np.arange(len(labels))
    fig_h = max(3.8, 0.32 * len(labels) + 1.2)
    fig, ax = plt.subplots(figsize=(6.2, fig_h), dpi=120)
    ax.scatter(-np.log10(pvals[~sig]), y[~sig], s=sizes[~sig], alpha=0.55, c="#64748b")
    ax.scatter(-np.log10(pvals[sig]), y[sig], s=sizes[sig], alpha=0.85, c="#0f766e")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("-log10(p-value)")
    ax.set_title("Pathway enrichment")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    log.info("pathway bubble %s", out.name)
    return str(out)


def integration_diagnostics(adata, batch_key: str, *, figdir: str | Path | None = None, max_cells: int = 8000) -> list[str]:
    """PCA/UMAP colored by batch, before vs after correction. Does not overwrite X_umap."""
    import numpy as np

    if batch_key not in adata.obs or int(adata.obs[batch_key].nunique()) < 2:
        return []
    out_dir = _figdir(figdir)
    idx = _subsample_idx(int(adata.n_obs), max_cells)
    batch = adata.obs[batch_key].to_numpy()[idx]
    saved: list[str] = []

    def _take(key: str):
        if key not in adata.obsm:
            return None
        x = np.asarray(adata.obsm[key])
        if x.ndim != 2 or x.shape[1] < 2:
            return None
        return x[idx, :2]

    pca_xy = _take("X_pca")
    p = _scatter_batch(pca_xy, batch, out_dir / "batch_pca_before.png", "PCA before correction (batch)") if pca_xy is not None else None
    if p:
        saved.append(str(p))
    after_key = next((k for k in ("X_pca_harmony", "X_scVI", "X_scanorama") if k in adata.obsm), None)
    if after_key:
        p = _scatter_batch(
            _take(after_key),
            batch,
            out_dir / "batch_pca_after.png",
            f"{after_key} after correction (batch)",
        )
        if p:
            saved.append(str(p))
    umap_xy = _take("X_umap")
    p = _scatter_batch(umap_xy, batch, out_dir / "batch_umap_after.png", "UMAP after correction (batch)") if umap_xy is not None else None
    if p:
        saved.append(str(p))
    if pca_xy is not None:
        try:
            from umap import UMAP

            pcs = np.asarray(adata.obsm["X_pca"])[idx]
            n_comp = min(30, pcs.shape[1])
            n_nb = min(15, max(2, len(idx) - 1))
            y = UMAP(n_components=2, n_neighbors=n_nb, min_dist=0.3, random_state=0).fit_transform(pcs[:, :n_comp])
            p = _scatter_batch(y, batch, out_dir / "batch_umap_before.png", "UMAP before correction (PCA→UMAP, batch)")
            if p:
                saved.append(str(p))
        except Exception as exc:
            log.info("uncorrected UMAP diagnostic skipped: %s", exc)
    return saved
