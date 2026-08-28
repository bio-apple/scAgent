# Yates et al. 2025 — Reconsidering mitochondrial QC cutoffs

Genome Biology. Challenges filtering all high-mito cells as dying.

## Key points

- High mito fraction may be technical (membrane rupture) or biological (metabolic tissue, some tumors, stress).
- Fixed 10% or 5% cutoffs systematically drop cells in tumor and metabolic studies.
- Correct approach: scatter `pct_counts_mt` vs `n_counts`; MAD outliers; tissue priors; report counts removed and rationale.

## Relation to OSCA / emptyDrops

- Cell calling: barcode rank / emptyDrops—not `n_genes>200` alone.
- QC trio: **Violin, Scatter, MAD decision** (required qc_expert outputs in this project).
