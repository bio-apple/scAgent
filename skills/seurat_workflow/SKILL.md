---
name: seurat-workflow
description: "scRNA core workflow: load matrix → QC → normalize → HVG → integrate → PCA → UMAP → Leiden. Seurat (R-first) or Scanpy backends."
capability: true
replaces: dataset_loader, qc_preprocessing, integration_batch, clustering_embedding, visualization, report_generation
---

# Capability: scRNA core workflow (Seurat / Scanpy)

## Scope
One capability covering **data intake through clustering** — not cell-type naming or confirmatory DE.

| Step | Python (Scanpy) | R-first |
|------|------------------|---------|
| Load | `read_h5ad` / `read_single_cell` (10x, loom, RDS) | Seurat + zellkonverter |
| QC | MAD/percentile per-sample; Scrublet | `pipeline_qc.R` MAD |
| Normalize | log1p / pearson / SCTransform (R fallback) | NormalizeData / SCTransform |
| HVG + PCA | seurat_v3 HVG, scale, PCA | FindVariableFeatures, ScaleData, RunPCA |
| Integrate | Harmony / scVI / BBKNN / Scanorama | Harmony (R) |
| Cluster | Leiden + joint resolution (silhouette + marker spread + size prior) | FindClusters |
| Viz | UMAP overview (leiden, batch, mito) | DimPlot |

## Honesty gates
- **Ambient**: auto → `none` until real SoupX wired; explicit `soupx` unavailable → counts unchanged
- **Doublet both**: second method may be `count_simulation` (labeled in metrics)
- **No fixed mito%<5**; per-sample MAD when `n_samples>1`
- **Integration skipped** when batch≡condition collinear

## Outputs
- `adata_qc.h5ad`, `qc_metrics.json`
- `.cache/after_cluster.h5ad`, `cluster_metrics.json` (split mode)
- `figures/violin_*`, `scatter_*`, `umap_overview.png`

## Scripts
- Phase 1: `workspace/qc_preprocess.py`
- Phase 2a (split): `workspace/cluster_only.py`
- Phase 2 (combined): `workspace/cluster_annotate.py`

## Related KB
- `knowledge/best_practices/qc.md`, `normalization.md`, `integration.md`, `clustering.md`
