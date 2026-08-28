# Heumos et al. 2023 — Best practices for single-cell analysis across modalities

Nature Reviews Genetics. Consensus workflow across single-cell modalities.

## Constraints for analysis agents

- Inspect distributions before QC thresholds. pctMT=10%, nHVG=2000 are not laws.
- Standard scRNA-seq: counts → QC (empty droplet / doublet / mito) → normalize → HVG → dim reduction → neighbors → cluster → annotate → DE.
- Group DE: biological replicates (samples) are units; cells are not independent. Exploratory markers: cell-level Wilcoxon; conclusions: pseudobulk + DESeq2/edgeR/limma.
- Integration is not default—needed for reference atlases; merge first when samples are already comparable.
- Spatial: declare resolution (spot / subcellular / single-cell). Spot-level claims ≠ single-cell facts.
- Report software versions, seeds, filter logs, reference dataset versions.

## When guidance breaks down

- Tiny or severely damaged libraries: unstable distributions; MAD may fail.
- Reference bias breaks every annotation layer (cross-species, developmental stage, tumor vs normal).
