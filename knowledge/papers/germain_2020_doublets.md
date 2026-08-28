# Germain et al. 2020 / Scrublet — doublet detection before annotation

Doublets create fake types (hybrid of two marker programs). Handle before annotation.

## Practice

- Run Scrublet / scDblFinder / SOLO in QC; scores in `obs`.
- Clusters with high `n_genes` + mixed markers: suspect doublets first.
- Do not annotate doublet clusters as novel cell types.
- High-loading 10x libraries: higher doublet rate; tune thresholds to expected rate.
