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
- FASTQ → counts (optional): nf-core/scrnaseq via archived wrapper recipes

## Outputs
- Primary object: **Seurat** (R-first) or **AnnData** (Python fallback)
- `obs`/`meta.data`: `sample`, optional `condition`, platform notes
- Detected **species** (human/mouse) and gene-id style (symbol vs Ensembl)
- Raw counts preserved (`layers['counts']` / `RNA@counts`)

## Procedure
1. Prefer Cell Ranger filtered outputs; do not assume 10x barcode rules for Parse Bio / other platforms.
2. Multi-path input → set `sample` per file; never silently merge without a batch key.
3. Keep raw counts; do not overwrite with log-normalized values.
4. Large h5ad: use backed mode when configured.
5. Reject processed-only matrices when downstream needs raw counts (DE, scVI) unless `layers['counts']` exists.

## Recipes

### R-first
```r
# Seurat::Read10X / Read10X_h5 → CreateSeuratObject
# Multi-sample: list of objects → merge(..., add.cell.ids = samples)
# h5ad ↔ Seurat: SeuratDisk / sceasy when needed
```

### Python
```python
import scanpy as sc
adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
# or sc.read_10x_h5 / sc.read_h5ad / sc.read_loom
adata.var_names_make_unique()
adata.layers["counts"] = adata.X.copy()
# multi-sample:
# adatas = [sc.read_10x_mtx(p) for p in paths]
# adata = sc.concat(adatas, label="sample", keys=sample_ids, index_unique="-")
```

### FASTQ (optional upstream)
- nf-core/scrnaseq samplesheet: `sample,fastq_1,fastq_2`
- Presets: `standard` (simpleaf) / `star` / `kallisto` / `cellranger`
- Velocity layers: STARsolo `--star-feature "Gene Velocyto"` or kb `lamanno|nac`
- Prefer combined filtered `.h5ad` as handoff into `qc_preprocessing`

## Gates
- Fail if counts are missing and only scaled matrix is present without a recovery path.
- Record n_cells, n_genes, n_samples in the run log.

## Related
- Next: `qc_preprocessing`
- Archive sources folded in: `nfcore-scrnaseq-wrapper`, `scrna-orchestrator` I/O rules
