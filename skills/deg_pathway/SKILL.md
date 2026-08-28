---
name: deg_pathway
description: "Cluster markers or condition DEG (pseudobulk when replicates≥2) plus GO/KEGG/GSEA interpretation."
---

# Scientific task: Differential expression & pathways

## Goal
Find genes and pathways that differ between clusters or biological conditions, with correct statistical units.

## DEG policy
| Design | Method |
|--------|--------|
| Cluster markers (exploratory) | Wilcoxon / MAST on cells — label as exploratory |
| Condition DE, **n_replicates ≥ 2** | **Pseudobulk + DESeq2/edgeR** (mandatory) |
| No biological replicates | Do not claim population inference; exploratory only |

## Pathway
- ORA / GSEA / GSVA (when libraries available)
- Prefer Hallmark / GO; report BH-FDR
- Gene-set choice matters more than the enrichment engine

## Outputs
- DEG table (gene, logFC, padj, direction)
- Volcano / heatmap of top genes
- Enrichment bar/dot plots
- Short biological interpretation (non-causal)

## Gates
- Never report raw p-values without multiple-testing correction for claims.
- Do not run cell-level Wilcoxon as the **group** conclusion when replicates exist.

## Related
- `knowledge/best_practices/pseudobulk-de.md`, `pathway-enrichment.md`, `marker-genes.md`

