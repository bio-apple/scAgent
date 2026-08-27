# Optional Monocle3 fate. scAgent does not require this package.
# CLI: Rscript monocle3.R input.h5ad status.json
# Writes status.json; never pretends a graph was learned if monocle3 is missing.

args <- commandArgs(trailingOnly = TRUE)
h5ad <- if (length(args) >= 1) args[[1]] else "workspace/adata_processed.h5ad"
out <- if (length(args) >= 2) args[[2]] else "workspace/monocle3_status.json"

.status <- function(status, detail = "") {
  txt <- sprintf('{"status":"%s","detail":"%s"}\n', status, gsub('"', "'", detail))
  writeLines(txt, out)
}

if (!requireNamespace("monocle3", quietly = TRUE)) {
  .status("skipped", "monocle3 not installed")
  quit(save = "no", status = 0)
}
if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  .status("skipped", "zellkonverter not installed; cannot read h5ad")
  quit(save = "no", status = 0)
}

sce <- zellkonverter::readH5AD(h5ad)
cds <- monocle3::new_cell_data_set(
  SingleCellExperiment::counts(sce),
  cell_metadata = as.data.frame(SummarizedExperiment::colData(sce)),
  gene_metadata = data.frame(gene_short_name = rownames(sce), row.names = rownames(sce))
)
cds <- monocle3::preprocess_cds(cds, num_dim = 20)
cds <- monocle3::reduce_dimension(cds)
cds <- monocle3::cluster_cells(cds)
cds <- monocle3::learn_graph(cds)
.status("ok", "learn_graph")
