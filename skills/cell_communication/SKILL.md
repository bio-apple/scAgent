---
name: cell_communication
description: "Ligand–receptor / CellChat or CellPhoneDB analysis on annotated cell types."
---

# Scientific task: Cell–cell communication

## Goal
Infer ligand–receptor interactions between **annotated** cell populations (hypothesis-generating).

## Prerequisites
- Stable cell type labels (`cell_annotation` first)
- Enough cells per type; filter tiny groups (e.g. min.cells ≥ 10)
- Expression should be **log-normalized** for CellPhoneDB-style tools (`X.max()` typically < ~10)

## Methods
| Tool | Role |
|------|------|
| **CellChat** (R) | Preferred in many publications; pathway + centrality |
| CellPhoneDB / LIANA | Python alternatives / consensus databases |
| NicheNet | Ligand → target gene regulatory focus (R) |
| MeboCost | Metabolite-mediated communication (optional specialty) |

## Recipes

### CellChat (R-first)
```r
library(CellChat)
data.input <- GetAssayData(obj, assay = "RNA", layer = "data")  # log-normalized
meta <- data.frame(labels = Idents(obj), row.names = colnames(obj))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat@DB <- CellChatDB.human  # or CellChatDB.mouse
cellchat <- subsetData(cellchat) |> identifyOverExpressedGenes() |> identifyOverExpressedInteractions()
cellchat <- computeCommunProb(cellchat, type = "triMean") |> filterCommunication(min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat) |> aggregateNet()
# netVisual_circle / bubble / chord; netAnalysis_computeCentrality
```

### CellPhoneDB via omicverse (Python)
```python
# Require categorical celltype_key without NA
cpdb_results, adata_cpdb = ov.single.cpdb_network_cal(
    adata, cpdb_file_path="cellphonedb.zip", celltype_key="cell_type",
    iterations=1000, threshold=0.1, pvalue=0.05)
# viz = ov.pl.CellChatViz(...); compute_aggregated_network; netVisual_circle / chord / bubble
```

### Metabolite (optional)
```python
# mebocost.create_obj(adata, group_col=..., species='human'|'mouse')
# mebo.infer_commu(n_permutations=1000); filter pval<0.05
```

## Outputs
- Interaction network / circle or hierarchy; L–R table with scores/p-values; sender–receiver notes

## Gates
- Hypothesis-generating only — not tissue proof of signaling.
- Prefer cell types over raw Leiden IDs when types exist.

## Related
- Folded from: CellChat cookbook, CellPhoneDB–omicverse mapping, metabolite-communication archive
