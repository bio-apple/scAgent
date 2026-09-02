---
name: visualization
description: "DEPRECATED → use seurat-workflow. UMAP/violin/dot plots."
deprecated: true
replaced_by: seurat-workflow
---

# Scientific task: Visualization

## Goal
Produce consistent, paper-ready figures from analysis objects — not a substitute for statistical tasks.

## Standard panels
- UMAP (sample, cluster, cell type, key genes)
- FeaturePlot / gene embedding
- Violin / VlnPlot
- DotPlot (markers / programs)
- Heatmap of markers or DE genes
- Optional: RidgePlot, PAGA overlay, communication circle/bubble (from upstream tasks)

## Defaults
- Clear fonts, minimal chrome (Nature-like)
- Export **PDF + PNG**, **300 dpi**
- Colorblind-friendly palettes when possible
- Before `color='leiden'|'cell_type'`: assert column exists (defensive)

## Recipes
```python
sc.pl.umap(adata, color=["sample", "leiden", "cell_type"], wspace=0.4)
sc.pl.dotplot(adata, var_names=markers, groupby="cell_type")
sc.pl.violin(adata, keys=["n_genes_by_counts", "pct_counts_mt"], groupby="sample")
# save: sc.settings.figdir; dpi=300; both .pdf and .png
```

```r
DimPlot(obj, group.by = c("sample", "seurat_clusters", "cell_type"))
FeaturePlot(obj, features = genes); DotPlot(obj, features = markers); VlnPlot(...)
ggsave(..., dpi = 300)
```

## Caveats (from DR literature)
- UMAP/t-SNE distances are not global metrics; do not claim hierarchy from 2D proximity alone.
- Prefer quantitative diagnostics (batch metrics, silhouette) over “looks mixed/separated”.

## Gates
- Every figure maps to an executed step (no fabricated panels).
- Titles/legends must match actual `obs` columns.

## Related
- Writer publication checklist; folded from orchestrator figure set + DR critique notes
