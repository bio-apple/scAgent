# Archived skills (not loaded at runtime)

Active product surface = **10 scientific-task skills** under `skills/*/SKILL.md`.

This folder keeps prior cookbooks for reference. **scRNA-relevant recipes were fused** into the active skills (2026-08):

| Archive (examples) | Fused into |
|--------------------|------------|
| nfcore-scrnaseq-wrapper | `dataset_loader` |
| single-cell-preprocessing-with-omicverse, scrna-orchestrator (QC parts) | `qc_preprocessing` |
| scrna-embedding, omicverse batch-correction | `integration_batch` |
| scrna-orchestrator / omicverse clustering | `clustering_embedding` |
| CellTypist / omicverse annotation / universal annotator | `cell_annotation` |
| orchestrator contrast DE | `deg_pathway` |
| Single-Cell Trajectory Inference | `trajectory` |
| CellChat + CellPhoneDB mapping + metabolite-communication | `cell_communication` |
| orchestrator figures / DR caveats | `visualization` |
| report contracts (orchestrator, embedding, nf-core) | `report_generation` |

## Still archived (out of active scRNA core)
- scATAC / Signac / ArchR
- Multiome / muon / CITE-seq WNN-heavy multimodal
- Perturb-seq / CRISPR screens
- Splicing / APA (non count-matrix core)
- Spatial / imaging / TCR-BCR packs (if present under legacy paths)

`legacy_granular/` may be absent on disk (git-deleted); recover with `git show HEAD:skills/_archive/...` when needed.
