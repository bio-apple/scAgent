---
name: trajectory
description: "Pseudotime / fate: PAGA, DPT, Palantir, scVelo (auto-gated). Monocle3 via R when requested."
capability: true
---

# Scientific task: Trajectory / fate

## Goal
Infer continuous structure or differentiation axes when biology is not purely discrete clusters.

## When / when not
- Use: differentiation, disease progression continua, ≥~200 cells, clusters already defined
- Skip: clearly discrete resting populations (e.g. resting PBMC panels) — do not force a fate axis

## Scope ladder
1. **Core**: PAGA + diffusion pseudotime (DPT) — needs root / early cell type
2. **+ Velocity**: scVelo only if `spliced`/`unspliced` layers exist (dynamical → fallback stochastic)
3. **+ Fate**: CellRank terminal probabilities (optional)
4. **R path**: Monocle3 / Slingshot when operating in Seurat

## Recipes

### Core (Python)
```python
sc.tl.paga(adata, groups="leiden")  # or cell type key
sc.pl.paga(adata, threshold=0.03)
# Set root: adata.uns['iroot'] = np.flatnonzero(adata.obs['cell_type'] == root)[0]
sc.tl.diffmap(adata)
sc.tl.dpt(adata)
```

### Velocity (only with layers)
```python
import scvelo as scv
scv.pp.moments(adata); scv.tl.velocity(adata, mode="dynamical")  # fallback: stochastic
scv.tl.velocity_graph(adata); scv.pl.velocity_embedding_stream(adata, basis="umap")
```

### Optional
- CellRank fate probabilities → `fate_probabilities.csv`
- StaVIA / VIA when user requests graph-based alternatives
- Gene–pseudotime trends + FDR for driver genes

## Outputs
- Pseudotime values; PAGA graph; gene-vs-pseudotime trends; velocity streams if run; confidence note (discrete vs continuous)

## Gates
- Do **not** force trajectory on discrete populations.
- Velocity without spliced/unspliced → skip and state why.
- Ask for input object + analysis scope (core / +velocity / +CellRank) before parameter deep-dives.

## Related
- Upstream: `clustering_embedding`; velocity only with spliced/unspliced layers
- Folded from: Single-Cell Trajectory Inference archive, StaVIA notes
