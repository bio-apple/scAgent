---
name: dataset_loader
description: "Load 10x/Cell Ranger, H5, H5AD, Loom, or Seurat RDS into AnnData/Seurat; detect species and gene symbols."
---

# Scientific task: Dataset loading

## Goal
Ingest raw or processed single-cell matrices and produce a single analysis object with consistent metadata.

## Inputs
- 10x `filtered_feature_bc_matrix/` (mtx/tsv/h5)
- `.h5` / `.h5ad` / `.loom`
- Seurat `.rds` / `.qs` (R path)
- Optional: sample sheet / condition labels

## Outputs
- Primary object: **Seurat** (R-first) or **AnnData** (Python fallback)
- `obs`/`meta.data`: `sample`, optional `condition`, platform notes
- Detected **species** (human/mouse) and gene-id style (symbol vs Ensembl)

## Procedure
1. Prefer Cell Ranger outputs; do not assume 10x barcode rules for Parse Bio / other platforms.
2. Multi-path input → set `sample` per file; never silently merge without a batch key.
3. Keep raw counts in a layer / assay (`counts` / `RNA@counts`); do not overwrite with log-normalized values.
4. Large h5ad: use backed mode when configured.

## R-first
```r
# Seurat::Read10X / Read10X_h5 → CreateSeuratObject
# SeuratDisk / sceasy for h5ad ↔ Seurat when needed
```

## Python fallback
```python
# scanpy.read_10x_mtx / read_10x_h5 / read_h5ad / read_loom
```

## Gates
- Fail if counts are missing and only scaled matrix is present without a recovery path.
- Record n_cells, n_genes, n_samples in the run log.

## Related
- Decision SOP: none (I/O). Upstream best practices apply after load → `qc_preprocessing`.

