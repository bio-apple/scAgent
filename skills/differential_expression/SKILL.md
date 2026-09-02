---
name: differential-expression
description: "Confirmatory pseudobulk DE (edgeR/DESeq2) and exploratory cluster markers (Wilcoxon). Never cell-level Wilcoxon as group conclusion."
capability: true
replaces: deg_pathway
---

# Capability: differential expression

## Two tiers (do not conflate)

| Tier | Question | Method | scAgent |
|------|----------|--------|---------|
| **Exploratory** | Which genes mark cluster X? | Wilcoxon / t-test (+ optional MAST) | `rank_genes` on `leiden` |
| **Confirmatory** | Group A vs B (with replicates)? | pseudobulk + edgeR/DESeq2 | `pseudobulk_de`; **hard-fail** if R unavailable |

## Confirmatory requirements
- ≥2 biological replicates per condition
- Aggregate: sample × cell_type (or cluster) raw counts
- Report FDR, effect size, engine used
- `confirmatory=True` without edgeR/DESeq2 → **RuntimeError** (not silent t-test)

## Exploratory cross-validation
- Optional second method (t-test vs Wilcoxon) → `cluster_marker_overlap.json`

## Outputs
- `cluster_marker_overlap.json`, rank_genes heatmap
- `pseudobulk_de` results in `uns` + metrics when condition comparison requested

## Script
- Split mode: `workspace/annotate_deg.py` (with annotation)
- Combined: tail of `cluster_annotate.py`

## Related
- `knowledge/best_practices/differential-expression.md`, `pseudobulk.md`
