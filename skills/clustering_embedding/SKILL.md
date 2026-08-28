---
name: clustering_embedding
description: "PCA → neighbors → Leiden/Louvain clusters → UMAP embedding and cluster stats."
---

# Scientific task: Clustering & embedding

## Goal
Build a neighborhood graph, define discrete clusters, and produce a 2D embedding for exploration.

## Pipeline (do not split into separate skills)
1. PCA (or use integrated latent if integration ran)
2. Neighbors / SNN graph
3. Leiden (or Louvain) clustering — scan resolutions; prefer silhouette / stability when available
4. UMAP (visualization only)

## Outputs
- Cluster labels in `obs` / `meta.data`
- PCA elbow / variance plot when useful
- UMAP colored by cluster and by sample/batch
- Cluster size table

## Gates
- Clustering uses **corrected** representation when integration was applied (`use_rep` / reduction).
- UMAP is not a statistical test; do not over-interpret distances.

## R-first
`RunPCA` → `FindNeighbors` → `FindClusters` → `RunUMAP`

## Python
`sc.tl.pca` → `sc.pp.neighbors` → `sc.tl.leiden` → `sc.tl.umap`

## Related
- `knowledge/best_practices/clustering.md`, `dimensionality-reduction.md`

