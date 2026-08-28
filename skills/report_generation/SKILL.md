---
name: report_generation
description: "Assemble Markdown/HTML report: summary, methods, QC numbers, limitations, session info."
---

# Scientific task: Report generation

## Goal
Turn pipeline artifacts into a reproducible analysis report suitable for collaborator review or supplementary methods.

## Must include
- Dataset & design summary (species, platform, n_samples, batch policy)
- QC decisions and filter rates
- Integration / clustering / annotation summary
- DEG & pathway highlights (if run)
- Literature-based best-practice notes (from RAG) when available
- Limitations (replicates, UMAP caveats, exploratory vs confirmatory)
- Session info / software versions / seed

## Formats
- Markdown (primary)
- HTML viewer when available
- Link notebooks / dual code–result blocks

## Gates
- Never describe plots or stats that were not generated.
- Label exploratory Wilcoxon vs pseudobulk confirmatory DE clearly.

## Related
- Writer agent (`agents/writer.py`); evidence chain section

