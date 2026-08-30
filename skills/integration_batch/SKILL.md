---
name: integration_batch
description: "Multi-sample batch integration (Harmony default; scVI/CCA when warranted). Skip if sample≡condition."
---

# Scientific task: Batch integration

## Goal
Correct technical batch effects across samples while preserving biology.

## When to run
- ≥2 samples / batches **and** sample is not 1:1 with biological condition.
- Skip when confounded (treatment ≡ batch) — prefer pseudobulk DE later, not overcorrection.

## Methods (priority)
| Situation | Method |
|-----------|--------|
| Default multi-sample | **Harmony** (on PCA) |
| Large (≥100k cells or ≥8 samples) | **scVI** (raw counts required) |
| Labels available for transfer | scANVI after scVI |
| R / Seurat workflows | Harmony or Seurat CCA / RPCA |
| Graph-level only | BBKNN (does not correct PCA) |

## Recipes

### Harmony (Python)
```python
sc.tl.pca(adata, n_comps=50)
sc.external.pp.harmony_integrate(adata, key="sample")  # or batch
sc.pp.neighbors(adata, use_rep="X_pca_harmony")
# then clustering_embedding uses this graph
```

### Harmony (R)
```r
obj <- RunPCA(obj) |> RunHarmony("sample") |> FindNeighbors(reduction = "harmony")
```

### scVI / scANVI (large / complex batches)
```python
# REQUIRE raw counts in layers['counts']; reject already-normalized-only input
# scvi.model.SCVI → obsm['X_scvi']; optional SCANVI with labels_key
# Downstream neighbors/UMAP: use_rep='X_scvi'
# Export integrated.h5ad: X_scvi + layers['counts'] preserved
```

### Diagnostics (not UMAP-alone)
- Prefer iLISI / kBET / ASW-batch / PCA batch R² when available
- Optional omicverse: `ov.single.batch_correction(..., methods='harmony'|'scVI')` + scib Benchmarker

## Outputs
- Corrected embedding (`X_pca_harmony` / `X_scvi` / integrated assay)
- Documented `sample_key` + method choice rationale

## Gates
- Forbidden: claim “integration succeeded” from UMAP mixing alone.
- scVI path must fail loudly if raw counts are missing.

## Related
- `knowledge/best_practices/integration.md`
- Folded from: `scrna-embedding`, Harmony cookbook, omicverse batch-correction skill
