# Seurat clustering + Azimuth annotation → h5ad with labels in obs
# Usage: Rscript pipeline_annotate.R <input_qc_h5ad> <output_h5ad> [tissue] [resolution]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript pipeline_annotate.R input.h5ad output.h5ad [tissue] [resolution]")
inp <- args[[1]]
out <- args[[2]]
tissue <- if (length(args) >= 3) args[[3]] else "default"
res <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.6

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed")
  if (!requireNamespace("zellkonverter", quietly = TRUE)) stop("zellkonverter not installed")
})

sce <- zellkonverter::readH5AD(inp)
obj <- Seurat::CreateSeuratObject(counts = SummarizedExperiment::assay(sce, 1))
if ("sample" %in% colnames(sce@colData)) obj$sample <- sce@colData$sample

obj <- Seurat::NormalizeData(obj, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
obj <- Seurat::ScaleData(obj, verbose = FALSE)
obj <- Seurat::RunPCA(obj, verbose = FALSE)
obj <- Seurat::FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- Seurat::FindClusters(obj, resolution = res, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, dims = 1:30, verbose = FALSE)

annotation_method <- "seurat_clusters"
obj$scagent_annotation <- as.character(obj$seurat_clusters)
obj$scagent_annotation_conf <- 0.5

if (requireNamespace("Azimuth", quietly = TRUE)) {
  tryCatch({
    az <- Azimuth::RunAzimuth(obj, reference = tissue)
    if ("predicted.celltype.l2" %in% colnames(az@meta.data)) {
      obj$scagent_annotation <- az$predicted.celltype.l2
      obj$scagent_annotation_conf <- az$predicted.celltype.l2.score
      annotation_method <- "azimuth"
    }
  }, error = function(e) message("Azimuth skipped: ", conditionMessage(e)))
}

if (identical(annotation_method, "seurat_clusters") && requireNamespace("SingleR", quietly = TRUE)) {
  tryCatch({
    sce_ann <- Seurat::as.SingleCellExperiment(obj)
    ref <- NULL
    if (requireNamespace("celldex", quietly = TRUE)) {
      ref <- celldex::HumanPrimaryCellAtlasData()
    }
    if (!is.null(ref)) {
      pred <- SingleR::SingleR(test = sce_ann, ref = ref, labels = ref$label.main)
      obj$scagent_annotation <- pred$labels
      obj$scagent_annotation_conf <- apply(pred$scores, 1, max)
      annotation_method <- "singler"
    }
  }, error = function(e) message("SingleR skipped: ", conditionMessage(e)))
}

sce_out <- Seurat::as.SingleCellExperiment(obj)
zellkonverter::writeH5AD(sce_out, out)

metrics <- list(
  engine = annotation_method,
  phase = "downstream",
  resolution = res,
  n_clusters = length(unique(obj$seurat_clusters)),
  n_cells = ncol(obj),
  scagent_annotation_method = annotation_method
)
write(jsonlite::toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), file = "r_metrics.json")
cat("SCAGENT_R_OK annotate method=", annotation_method, "\n", sep = "")
