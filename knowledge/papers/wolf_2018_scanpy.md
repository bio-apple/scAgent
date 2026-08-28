# Wolf et al. 2018 — SCANPY: large-scale single-cell gene expression data analysis

Genome Biology. Scalable Python scRNA-seq framework centered on AnnData.

## Skills mapping in this repo

- Standard exploration: `skills/scanpy-scrna-seq`
- Data structures: `skills/anndata-data-structure`
- Default clustering: Leiden (Traag 2019)—do not use Louvain without reason.
- Neighbor graph → UMAP is visualization; cluster on the graph, not UMAP coordinates.

## Platform notes

- 10x Genomics: `sc.read_10x_mtx` / `sc.read_10x_h5`, `var_names='gene_symbols'`.
- Parse Biosciences: different barcode structure—planner marks `platform=parse`, not 10x chemistry assumptions.
- Multi-sample: concat with `sample` column preserved, then decide Harmony/scVI.
