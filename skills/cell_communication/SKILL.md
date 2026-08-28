---
name: cell_communication
description: "Ligand–receptor / CellChat or CellPhoneDB analysis on annotated cell types."
---

# Scientific task: Cell–cell communication

## Goal
Infer ligand–receptor interactions between annotated cell populations.

## Prerequisites
- Stable **cell type** labels (run `cell_annotation` first)
- Sufficient cells per type; filter tiny clusters

## Methods
- **CellChat** (preferred in R / many publications)
- CellPhoneDB / LIANA as alternatives

## Outputs
- Interaction network / circle or hierarchy plot
- L–R pair table with scores / p-values
- Notes on sender–receiver populations of interest

## Gates
- Communication is **hypothesis-generating**, not proof of signaling in tissue.
- Do not run on unannotated Leiden IDs alone when types are available.

## Related
- Literature RAG / CellChat skill archive for R details

