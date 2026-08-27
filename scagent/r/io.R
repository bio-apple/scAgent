# Convert Seurat RDS <-> h5ad via zellkonverter (optional R dependency).
# Usage: Rscript io.R read input.rds output.h5ad
#        Rscript io.R write input.h5ad output.rds

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: Rscript io.R read|write input output")
}
cmd <- args[[1]]
inp <- args[[2]]
out <- args[[3]]

if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  stop("Install R package zellkonverter (and Seurat for RDS input)")
}

if (cmd == "read") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Install Seurat to read .rds")
  }
  obj <- readRDS(inp)
  if (inherits(obj, "Seurat")) {
    sce <- Seurat::as.SingleCellExperiment(obj)
  } else {
    sce <- obj
  }
  zellkonverter::writeH5AD(sce, out)
} else if (cmd == "write") {
  sce <- zellkonverter::readH5AD(inp)
  saveRDS(sce, out)
} else {
  stop("unknown command: ", cmd)
}
