---
name: spatial-analysis
description: "Spatial transcriptomics (Visium / Xenium / ST): deconvolution, domain detection, spatial viz. Placeholder — not on default scRNA path."
capability: true
status: placeholder
---

# Capability: spatial analysis (ST)

## Status: **placeholder**
scAgent default pipeline is **scRNA matrix → AnnData**. Spatial modules are archived under `skills/_archive/` until a Squidpy/Giotto path is wired into templates.

## Intended scope (future)
| Task | Tools |
|------|-------|
| Load | Squidpy, Giotto, Seurat spatial |
| QC | spot-level counts, mitochondrial %, spatial outliers |
| Deconvolution | cell2location, RCTD, stereoscope |
| Domains | Banksy, SpatialDE, BayesSpace |
| Viz | spatial scatter, co-occurrence, neighborhood enrichment |

## When user asks for ST
1. Confirm platform (Visium HD, Xenium, Slide-seq, etc.)
2. Warn: reference scRNA atlas must **match tissue** (resolution honesty)
3. Route to archived cookbooks or external nf-core spatial until implemented

## Related archive skills
- `skills/_archive/` spatial-deconvolution, spatial-visualization, spatial-domains
