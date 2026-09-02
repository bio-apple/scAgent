---
name: clustering_embedding
description: "DEPRECATED → use seurat-workflow. PCA, Leiden, UMAP."
deprecated: true
replaced_by: seurat-workflow
---

# Scientific task: Clustering & embedding

## Goal
Build a neighborhood graph, define discrete clusters, and produce a 2D embedding for exploration.

## Pipeline (do not split into separate skills)
1. PCA — **or** use integrated latent if `integration_batch` ran (`X_pca_harmony` / `X_scvi`)
2. Neighbors / SNN graph on the chosen representation
3. Leiden (or Louvain) — scan resolutions; prefer silhouette / stability when available
4. UMAP (visualization only)

## Recipes

### R-first
```r
# If Harmony: reduction = "harmony"; else "pca"
obj <- FindNeighbors(obj, reduction = red, dims = 1:30) |>
  FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8)) |>
  RunUMAP(reduction = red, dims = 1:30)
```

### Python
```python
use_rep = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else (
    "X_scvi" if "X_scvi" in adata.obsm else "X_pca"
)
if use_rep == "X_pca" and "X_pca" not in adata.obsm:
    sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, use_rep=use_rep)
sc.tl.leiden(adata, resolution=0.5)  # sweep 0.2–1.0 as needed
sc.tl.umap(adata)
```

### Optional alternatives (when requested)
- omicverse: `ov.utils.cluster(method='leiden'|'scICE'|'GMM')`; cNMF for programs
- GPU neighbors: `method='cagra'` only if RAPIDS stack is available
- Defensive: before plotting `color='leiden'`, ensure neighbors + leiden exist

## Outputs
- Cluster labels in `obs` / `meta.data`
- PCA elbow / variance plot when useful
- UMAP by cluster **and** sample/batch
- Cluster size table

## Gates
- Clustering uses **corrected** representation when integration was applied.
- UMAP is not a statistical test; do not over-interpret distances.
- Do not re-run QC/doublets here.

## Related
- `knowledge/best_practices/clustering.md`, `dimensionality-reduction.md`
- Folded from: `scrna-orchestrator`, omicverse clustering/batch skill
