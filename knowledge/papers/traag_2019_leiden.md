# Traag, Waltman & van Eck 2019 — From Louvain to Leiden

Scientific Reports. Leiden fixes Louvain disconnected-community issues.

## Clustering constraints

- Default Leiden. Resolution is not fixed at 0.8—calibrate with stability, marker separability, biological granularity.
- Over-clustering: one type split into many clusters → reviewer warns UMAP/Leiden too fine.
- Under-clustering: distinct lineages merged → raise resolution or hierarchical clustering.
- On developmental continua, discrete clusters are approximate; use PAGA / pseudotime and state assumptions.
