# Seurat QC + normalize → h5ad (scAgent R-first path)
# Usage: Rscript pipeline_qc.R <input> <output_h5ad> [sample_key] [nmads]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript pipeline_qc.R input output.h5ad [sample_key] [nmads]")
inp <- args[[1]]
out <- args[[2]]
sample_key <- if (length(args) >= 3) args[[3]] else "sample"
nmads <- if (length(args) >= 4) as.numeric(args[[4]]) else 5

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed")
  if (!requireNamespace("zellkonverter", quietly = TRUE)) stop("zellkonverter not installed")
})

n_before <- NA_integer_
n_after <- NA_integer_

load_seurat <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) stop("RDS is not a Seurat object")
    return(obj)
  }
  if (grepl("\\.h5ad$", path, ignore.case = TRUE)) {
    sce <- zellkonverter::readH5AD(path)
    return(Seurat::CreateSeuratObject(counts = SummarizedExperiment::assay(sce, 1)))
  }
  stop("input must be .rds or .h5ad")
}

obj <- load_seurat(inp)
n_before <- ncol(obj)
obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
obj <- subset(obj, subset = nFeature_RNA > 0)
med <- median(obj$nFeature_RNA)
mad <- mad(obj$nFeature_RNA)
upper <- med + nmads * mad
lower <- max(0, med - nmads * mad)
obj <- subset(obj, subset = nFeature_RNA > lower & nFeature_RNA < upper & percent.mt < median(obj$percent.mt) + nmads * mad(obj$percent.mt))
n_after <- ncol(obj)
obj <- Seurat::NormalizeData(obj, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- Seurat::ScaleData(obj, verbose = FALSE)
obj <- Seurat::RunPCA(obj, verbose = FALSE)
sce <- Seurat::as.SingleCellExperiment(obj)
zellkonverter::writeH5AD(sce, out)

metrics <- list(
  engine = "seurat",
  phase = "qc",
  n_before = n_before,
  n_after = n_after,
  pct_removed = if (!is.na(n_before) && n_before > 0) round(100 * (1 - n_after / n_before), 2) else 0,
  nmads = nmads,
  qc_method = "mad_seurat"
)
write(jsonlite::toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), file = "r_metrics.json")
cat("SCAGENT_R_OK qc n_before=", n_before, " n_after=", n_after, "\n", sep = "")
