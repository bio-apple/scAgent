# Korsunsky et al. 2019 — Harmony: fast, sensitive and accurate integration

Nature Methods. Iterative soft clustering + within-cluster linear correction in PCA space.

## When it fits

- Input: PCA embedding, not raw counts. Output: `obsm['X_pca_harmony']`.
- Does not modify expression matrix. Neighbors / UMAP / Leiden use corrected embedding.
- Fast to ~million cells; multiple covariates (batch, donor, platform).

## Limits

- Linear correction. Strong nonlinear batch or no shared types: overcorrect or fail.
- Cannot fix bad design: if treatment is fully collinear with batch, Harmony removes treatment as “batch”.

## vs scVI

- Fast, reproducible, simple batch → Harmony.
- Generative model, imputation, multimodal, uncertainty DE → scVI/scANVI (Lopez 2018; Gayoso 2022 scvi-tools).
