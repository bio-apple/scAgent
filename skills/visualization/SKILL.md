---
name: visualization
description: "Publication-style figures: UMAP, Feature/Violin/Dot/Heatmap; PDF+PNG @ 300 dpi."
---

# Scientific task: Visualization

## Goal
Produce consistent, paper-ready figures from analysis objects — not a substitute for statistical tasks.

## Standard panels
- UMAP (sample, cluster, cell type, key genes)
- FeaturePlot / FeaturePlot-like
- Violin / VlnPlot
- DotPlot
- Heatmap of markers or DE genes
- Optional RidgePlot

## Defaults
- Nature-like theme (clear fonts, minimal chrome)
- Export **PDF + PNG**, **300 dpi**
- Colorblind-friendly palettes when possible

## Gates
- Every figure must map to an executed analysis step (no fabricated panels).
- Titles/legends must match actual `obs` columns.

## Related
- Publication figure checklist in the Writer agent

