---
name: report_generation
description: "Assemble Markdown/HTML report: summary, methods, QC numbers, limitations, session info."
---

# Scientific task: Report generation

## Goal
Turn pipeline artifacts into a reproducible analysis report for collaborator review or methods supplements.

## Must include
- Dataset & design (species, platform, n_samples, batch policy)
- QC decisions and filter rates (per sample)
- Integration / clustering / annotation summary (+ method rationale)
- DEG & pathway highlights if run — label **exploratory** vs **pseudobulk confirmatory**
- Trajectory / communication summaries if run
- Literature-based best-practice notes (RAG) when available
- Limitations (replicates, UMAP caveats, confounded batches)
- Session info / software versions / seed / key file checksums when available

## Artifact contract (align with orchestrator-style outputs)
- `report.md` (primary) + optional HTML
- Tables: QC summary, cluster sizes, markers/DEG, optional L–R
- Figures linked with captions that match real files
- Optional `result.json` / provenance for machine handoff

## Formats
- Markdown primary; HTML viewer when available; link notebooks / dual code–result blocks

## Gates
- Never describe plots or stats that were not generated.
- Clear exploratory vs confirmatory language for DE.

## Related
- Writer agent (`agents/writer.py`); evidence chain
- Folded from: scrna-orchestrator / embedding / nf-core report contracts
