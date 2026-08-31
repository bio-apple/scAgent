---
name: deg_pathway
description: "Cluster markers or condition DEG (pseudobulk when replicates≥2) plus GO/KEGG/GSEA interpretation."
---

# Scientific task: Differential expression & pathways

## Goal
Find genes and pathways that differ between clusters or conditions, with correct statistical units.

## DEG policy
| Design | Method |
|--------|--------|
| Cluster markers (exploratory) | Wilcoxon / MAST on cells — label **exploratory** |
| Condition DE, **n_replicates ≥ 2** | **Pseudobulk + DESeq2/edgeR** (mandatory) |
| No biological replicates | Exploratory only — no population inference claims |
| Within-cluster contrasts | Optional all-pairs Wilcoxon inside a partition (exploratory) |

## Recipes

### Cluster markers (Python)
```python
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", pts=True)
sc.get.rank_genes_groups_df(adata, group=None)
```

### Pseudobulk condition DE (when replicates ≥ 2)
```python
# Aggregate counts per (sample × cell_type); DESeq2/edgeR in R (preferred) or pyDESeq2
# Never use cell-level Wilcoxon as the group-level conclusion when replicates exist
```

### R
```r
# Exploratory: FindAllMarkers / FindMarkers (Wilcoxon/MAST)
# Confirmatory: AggregateExpression / muscat → DESeq2 or edgeR on sample-level counts
```

### Pathway
- ORA / GSEA / GSVA; prefer Hallmark / GO; report BH-FDR
- Gene-set choice matters more than the enrichment engine

## Outputs
- DEG table (gene, logFC, padj, direction); volcano / heatmap; enrichment plots; short non-causal interpretation

## Gates
- Condition DE requires **annotated cell types** (or stable cluster→type map) before `pseudobulk_de`.
- `condition_key` + **n_replicates ≥ 2** → mandatory pseudobulk + DESeq2/edgeR; exploratory DEG intent alone is not confirmatory.
- Never invent a condition column (`unspecified`) to force DE to run.
- No raw p-values as claims without multiple-testing correction.
- Label exploratory vs confirmatory clearly in report text.

## Related
- `knowledge/best_practices/pseudobulk-de.md`, `pathway-enrichment.md`
- Folded from: scrna-orchestrator contrast modes, Seurat core DE notes
