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
| Large (≥100k cells or ≥8 samples) | **scVI** |
| R / Seurat workflows | Harmony or Seurat CCA / RPCA |
| Graph-level only | BBKNN (does not correct PCA) |

## Outputs
- Corrected embedding (`harmony` / `X_scVI` / integrated assay)
- Before/after UMAP is optional visualization — **not** proof of success
- Metrics when available: iLISI, kBET, PCA batch R²

## Gates
- Forbidden: claim “integration succeeded” from UMAP mixing alone.
- Document `sample_key` and why method was chosen.

## Related
- `knowledge/best_practices/integration.md`

