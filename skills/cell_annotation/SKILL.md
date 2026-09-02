---
name: cell-annotation
description: "Assign cell types via reference mapping + markers with dual validation. Requires --tissue for catalog/CellTypist."
capability: true
---

# Scientific task: Cell type annotation

## Goal
Label clusters/cells with cell types using layered evidence — never a single marker gene alone.

## Strategy (auto-select)
| Context | Primary tools |
|---------|----------------|
| PBMC / immune | **CellTypist** (`Immune_All_Low.pkl`) + majority vote; R: Azimuth / SingleR |
| Has reference atlas | SingleR / scANVI / popV / weighted kNN transfer |
| Unknown tissue | Cluster markers ∩ catalog (CellMarker / Panglao) + dual validation |
| Multi-source labels | Consensus vote (e.g. CellVote / majority across methods) |

## Required practice
- Dual validation: reference label **and** ≥2 positive / ≥1 negative markers
- Hierarchical labels when possible (lineage → subtype)
- Low-confidence → `unvalidated` / `mixed` — do not force a type
- Cite marker sources (`marker_db`, ontology) in the report

## Recipes

### CellTypist (Python, immune-first)
```python
import celltypist
from celltypist import models
# Prefer local/cached models in production (avoid runtime download when offline)
pred = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True)
adata = pred.to_adata()  # predicted_labels, majority_voting, conf_score
# Flag conf_score < 0.5 for manual review
```

### Marker + catalog (always as validation)
```python
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", pts=True)
# Map top markers to knowledge/marker_db + tissue panels; require multi-gene support
```

### Optional (omicverse)
```python
# SCSA: clustertype='leiden' (NOT cluster=); target='cellmarker'|'panglaodb'
# scsa = ov.single.pySCSA(adata, foldchange=1.5, pvalue=0.01, species='Human', tissue='All', target='cellmarker')
# scsa.cell_anno(clustertype='leiden'); scsa.cell_auto_anno(adata, clustertype='leiden', key='scsa_celltype')
# COSG writes adata.uns['rank_genes_groups'] — does NOT create obs celltype columns
# CellVote / GPTAnno: only with explicit user request + credentials; still require marker validation
```

### R
```r
# SingleR::SingleR(test, ref, labels=ref$label.main)
# Azimuth for supported tissues; Azimuth/SingleR labels still need marker DotPlot check
```

## Outputs
- Cell type UMAP; marker DotPlot; confidence / agreement summary

## Gates
- Dual validation is **expression-gated** (`dual_validate_expression`): ≥2 positive genes above threshold and ≥1 negative below — not catalog list length alone.
- R/Azimuth/SingleR success is **one fuse source only**; never `SystemExit` past marker dual + `fuse_annotation`.
- Unknown / empty tissue: refuse silent PBMC marker fallback; mouse: no human CellTypist/catalog without homology.
- Cluster IDs (`leiden` / `seurat_clusters`) are never final cell types; write ontology IDs when catalog has `cl_id`.

## Related
- `knowledge/best_practices/cell-annotation.md`, `marker-genes.md`; KB: `marker_db`, `cell_ontology`
- Folded from: CellTypist cookbook, omicverse annotation, universal annotator, orchestrator CellTypist path
