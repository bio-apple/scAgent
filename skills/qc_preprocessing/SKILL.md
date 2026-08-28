---
name: qc_preprocessing
description: "Tissue-aware QC, doublet/ambient handling, normalize, HVG, scale — required before clustering."
---

# Scientific task: QC & preprocessing

## Goal
Remove empty droplets / dying cells / doublets, normalize, select HVGs, and prepare a scaled matrix for PCA — **without** entering clustering.

## Must include
1. QC metrics: n_genes, total counts, mito % (jointly; per-sample)
2. Filter with **MAD / percentile** (not fixed mito%<5 by default)
3. Doublet detection (Scrublet ± scDblFinder); optional ambient (SoupX/CellBender)
4. Normalize → log1p → HVG → Scale (or equivalent Seurat pipeline)

## Required figures / tables
- Violin (QC metrics)
- Scatter (counts vs genes / mito)
- QC summary table (n kept/removed per sample)

## Gates
- **Do not proceed** to integration/clustering if QC plots/summary are missing.
- Warn if >~30% cells removed (configurable).
- Never default to hard mito%<5 across tissues (heart/kidney/tumor need profiles).

## R-first
`PercentageFeatureSet` → subset → `NormalizeData` → `FindVariableFeatures` → `ScaleData`

## Python
`sc.pp.calculate_qc_metrics` → MAD filter → `normalize_total`/`log1p` → `highly_variable_genes` → `scale`

## Related
- `knowledge/best_practices/qc.md`, `doublet-detection.md`, `normalization.md`, `feature-selection.md`

