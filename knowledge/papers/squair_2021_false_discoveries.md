# Squair et al. 2021 — Confronting false discoveries in single-cell differential expression

Nature Communications. Cell-level DE treats cells as independent replicates → severe false positives.

## Core rules

- Biological replicates are samples. Five thousand T cells from one donor are not 5,000 independent observations.
- Group conclusions (disease vs control, treatment vs control): **pseudobulk**—sum counts per sample × cell type, then DESeq2 / edgeR / limma-voom.
- Cell-level Wilcoxon / t-test: exploratory markers and within-cluster enrichment only—not primary paper p-values.
- Crowell et al. 2020 (muscat): multi-sample, multi-cell-type pseudobulk framework.

## Agent review checklist

- Does code aggregate by `sample`?
- Is replicate count N reported?
- Multiple testing correction (FDR / BH)?
- Single-sample studies must state “exploratory; no between-group inference” in the report.
