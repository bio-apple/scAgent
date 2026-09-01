---
name: qc_preprocessing
description: "Tissue-aware QC, doublet/ambient handling, normalize, HVG, scale — required before clustering."
---

# Scientific task: QC & preprocessing

## Goal
Remove empty droplets / dying cells / doublets, normalize, select HVGs, and prepare a scaled matrix for PCA — **without** entering clustering.

## Must include
1. QC metrics: n_genes, total counts, mito % (jointly; **per-sample**)
2. Filter with **MAD / percentile** (not fixed mito%<5 by default)
3. Doublet detection (Scrublet / scDblFinder / DoubletFinder); ambient only via real SoupX/DecontX (heuristic must be explicit)
4. Normalize → log1p → HVG → Scale (or SCTransform / Pearson residuals)

## Decision heuristics
| Tissue / context | Mito starting point | Notes |
|------------------|---------------------|-------|
| PBMC / blood | ~5–10% | Also track HB genes |
| Tumor / stressed | ~15–20% | MAD preferred over hard cut |
| Heart / kidney | higher mito OK | Tissue-aware profiles |
| Doublets | run **per sample** on raw counts | Scrublet before normalize |

- Ribosomal (RPS/RPL) and hemoglobin optional QC covariates.
- If >~30% cells removed → warn and re-check thresholds.
- HVG: prefer `seurat_v3` on counts; if LOESS fails (small batches <500) → fallback `seurat` / `cell_ranger`.
- Always `store` raw counts before transforms (`layers['counts']`).

## Recipes

### R-first
```r
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")  # mouse: "^mt-"
# MAD or percentile gates per sample — avoid global mito<5
# DoubletFinder / scDblFinder per sample, then subset
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData()
# Alternative: SCTransform(obj, vars.to.regress = "percent.mt")
```

### Python (scanpy)
```python
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
# MAD filter (nmads≈3–5) per sample, then:
# scrublet / sc.pp.scrublet on counts BEFORE normalize
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
sc.pp.scale(adata, max_value=10)
```

### Optional (omicverse)
```python
# ov.pp.qc(..., tresh={mito_perc, nUMIs, detected_genes}, doublets_method='scrublet')
# ov.utils.store_layers(adata, layers='counts')
# ov.pp.preprocess(adata, mode='shiftlog|pearson', n_HVGs=2000)
# batch col: fillna → category before any batch-aware step
```

## Required figures / tables
- Violin (QC metrics), scatter (counts vs genes / mito)
- QC summary table (n kept/removed per sample)

## Gates
- **Do not proceed** to integration/clustering if QC plots/summary are missing.
- Never default to hard mito%<5 across tissues.
- Scrublet/doublet calls belong here — not inside clustering skill.

## Related
- `knowledge/best_practices/doublet-detection.md`, `normalization.md`, `feature-selection.md`
- Folded from: `single-cell-preprocessing-with-omicverse`, `scrna-orchestrator`, Seurat core QC, doublet guides
