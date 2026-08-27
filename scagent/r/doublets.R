# scDblFinder on an h5ad; write barcode,score,class CSV.
# Usage: Rscript doublets.R input.h5ad output.csv [sample_column|NONE]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript doublets.R input.h5ad output.csv [sample_col|NONE]")
}
inp <- args[[1]]
out <- args[[2]]
sample_col <- if (length(args) >= 3) args[[3]] else "NONE"

if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  stop("Install R package scDblFinder")
}
if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  stop("Install R package zellkonverter")
}

sce <- zellkonverter::readH5AD(inp)
if (sample_col != "NONE" && sample_col %in% colnames(SummarizedExperiment::colData(sce))) {
  sce <- scDblFinder::scDblFinder(sce, samples = sample_col)
} else {
  sce <- scDblFinder::scDblFinder(sce)
}
cd <- SummarizedExperiment::colData(sce)
cls <- as.character(cd$scDblFinder.class)
score <- as.numeric(cd$scDblFinder.score)
df <- data.frame(
  barcode = colnames(sce),
  score = score,
  class = cls,
  stringsAsFactors = FALSE
)
utils::write.csv(df, out, row.names = FALSE)
