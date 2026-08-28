# Normalization: LogNormalize vs scran vs SCTransform

- **LogNormalize** (Seurat / Scanpy `normalize_total` + `log1p`): simple, reproducible; enough for most exploratory work.
- **scran deconvolution**: pooled size factors; robust to composition shifts (Lun / Bioconductor).
- **SCTransform**: regularized NB regression for technical noise; narrative differs, often correlates with log results.

Default Python path: `normalize_total(target_sum=1e4)` + `log1p` (CellTypist convention). Switch to scran/SCT only when the user asks.
