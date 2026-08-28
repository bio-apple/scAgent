# Integration is a decision, not a default

| Scenario | Recommendation |
|----------|----------------|
| Multi-sample technical batch, shared cell types | Harmony, then scIB dual-metric check |
| Nonlinear batch / generative model needed | scVI |
| Homogeneous biological replicates, DE-focused | Merge first; batch as covariate |
| Treatment fully collinear with batch | Do not “remove batch”—you remove treatment signal |

UMAP mixing ≠ success. Overcorrection erases disease differences.
