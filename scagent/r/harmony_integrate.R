# Harmony integration on Seurat object (R)
# Usage: Rscript harmony_integrate.R <input_h5ad> <output_h5ad> <batch_key>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript harmony_integrate.R input.h5ad output.h5ad batch_key")
inp <- args[[1]]
out <- args[[2]]
batch_key <- args[[3]]

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed")
  if (!requireNamespace("harmony", quietly = TRUE)) stop("harmony R package not installed")
  if (!requireNamespace("zellkonverter", quietly = TRUE)) stop("zellkonverter not installed")
})

sce <- zellkonverter::readH5AD(inp)
obj <- Seurat::CreateSeuratObject(counts = SummarizedExperiment::assay(sce, 1))
for (col in colnames(sce@colData)) obj[[col]] <- sce@colData[[col]]

obj <- Seurat::NormalizeData(obj, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
obj <- Seurat::ScaleData(obj, verbose = FALSE)
obj <- Seurat::RunPCA(obj, verbose = FALSE)
if (!(batch_key %in% colnames(obj@meta.data))) stop("batch_key not in metadata: ", batch_key)
obj <- harmony::RunHarmony(obj, group.by.vars = batch_key, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)

sce_out <- Seurat::as.SingleCellExperiment(obj)
zellkonverter::writeH5AD(sce_out, out)

metrics <- list(engine = "harmony_r", phase = "integration", integrator = "harmony_r")
write(jsonlite::toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), file = "r_metrics.json")
cat("SCAGENT_R_OK harmony\n")
