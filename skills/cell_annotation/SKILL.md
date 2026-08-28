---
name: cell_annotation
description: "Assign cell types via reference mapping + markers with dual validation (not single-gene labels)."
---

# Scientific task: Cell type annotation

## Goal
Label clusters/cells with cell types using layered evidence.

## Strategy (auto-select)
| Context | Primary tools |
|---------|----------------|
| PBMC / immune | CellTypist (+ Azimuth in R) |
| Has reference atlas | SingleR / scANVI / popV |
| Unknown tissue | Marker panels + cluster DE ∩ catalog |

## Required practice
- Dual validation: reference label **and** ≥2 positive / ≥1 negative markers
- Hierarchical labels when possible (lineage → subtype)
- Low-confidence cells marked `unvalidated` / `mixed` — do not force a type

## Outputs
- Cell type UMAP
- Marker DotPlot / DotPlot of top markers
- Confidence / agreement summary

## Gates
- Forbid single-gene annotation as final label.
- Cite marker sources (catalog / ontology) in the report.

## Related
- `knowledge/best_practices/cell-annotation.md`, `marker-genes.md`
- Structured KB: `marker_db`, `cell_ontology`, `tissue_reference`

