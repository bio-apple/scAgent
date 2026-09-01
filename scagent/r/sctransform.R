# Optional Seurat SCTransform → h5ad for scAgent.
# Usage: Rscript sctransform.R <input.h5ad> <output.h5ad>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript sctransform.R input.h5ad output.h5ad")
inp <- args[[1]]
out <- args[[2]]

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed")
  if (!requireNamespace("zellkonverter", quietly = TRUE)) stop("zellkonverter not installed")
})

sce <- zellkonverter::readH5AD(inp)
obj <- Seurat::CreateSeuratObject(counts = SummarizedExperiment::assay(sce, 1))
obj <- Seurat::SCTransform(obj, verbose = FALSE)
# Write SCT assay as primary assay for AnnData.X (Pearson-like residuals)
DefaultAssay(obj) <- "SCT"
sce_out <- Seurat::as.SingleCellExperiment(obj, assay = "SCT")
zellkonverter::writeH5AD(sce_out, out)
cat("SCAGENT_R_OK sctransform n_cells=", ncol(obj), "\n", sep = "")
