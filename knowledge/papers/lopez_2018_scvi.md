# Lopez et al. 2018 — Deep generative modeling for single-cell transcriptomics (scVI)

Nature Methods. VAE on UMI counts (ZINB/NB); latent space separates biology and batch.

## Practice

- Input must be **raw counts**, not log-normalized. Store in `layers['counts']`.
- `setup_anndata` → `SCVI` → `train` → `get_latent_representation`.
- scANVI: semi-supervised label transfer; totalVI: CITE-seq; DestVI: spatial deconvolution.
- Post-2025 foundation-model benchmarks: simple baselines often match scFM. Validate new models on user data first.

## When not to use

- <10k cells, simple batch, need minutes: Harmony is enough.
- No GPU and very large n: training cost is high.
