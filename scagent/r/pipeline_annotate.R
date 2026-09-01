# Seurat clustering + Azimuth/SingleR annotation → h5ad with labels in obs
# Usage: Rscript pipeline_annotate.R <input_qc_h5ad> <output_h5ad> [tissue] [resolution]
# NOTE: Reference labels alone are NOT final cell types — Python evidence layer must fuse.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript pipeline_annotate.R input.h5ad output.h5ad [tissue] [resolution]")
inp <- args[[1]]
out <- args[[2]]
tissue <- if (length(args) >= 3) tolower(args[[3]]) else "default"
res <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.6

# Map scAgent tissue keys → Azimuth reference names (SeuratData / Azimuth hubs).
azimuth_ref_map <- c(
  pbmc = "pbmcref",
  blood = "pbmcref",
  immune = "pbmcref",
  lung = "lungref",
  airway = "lungref",
  bone = "bonemarrowref",
  marrow = "bonemarrowref",
  brain = "humancortexref",
  cortex = "humancortexref"
)
azimuth_ref <- unname(azimuth_ref_map[tissue])
if (length(azimuth_ref) == 0 || is.na(azimuth_ref)) azimuth_ref <- NULL

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

# Never treat cluster IDs as cell types — leave unassigned until a real mapper succeeds.
annotation_method <- "unassigned"
obj$scagent_annotation <- "unassigned"
obj$scagent_annotation_conf <- 0.0

if (!is.null(azimuth_ref) && requireNamespace("Azimuth", quietly = TRUE)) {
  tryCatch({
    az <- Azimuth::RunAzimuth(obj, reference = azimuth_ref)
    if ("predicted.celltype.l2" %in% colnames(az@meta.data)) {
      obj$scagent_annotation <- az$predicted.celltype.l2
      obj$scagent_annotation_conf <- az$predicted.celltype.l2.score
      annotation_method <- "azimuth"
    } else if ("predicted.celltype.l1" %in% colnames(az@meta.data)) {
      obj$scagent_annotation <- az$predicted.celltype.l1
      obj$scagent_annotation_conf <- az$predicted.celltype.l1.score
      annotation_method <- "azimuth"
    }
  }, error = function(e) message("Azimuth skipped: ", conditionMessage(e)))
} else if (is.null(azimuth_ref)) {
  message("Azimuth skipped: no mapped reference for tissue=", tissue)
}

if (identical(annotation_method, "unassigned") && requireNamespace("SingleR", quietly = TRUE)) {
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
# Preserve leiden-compatible cluster column for Python fusion
sce_out$leiden <- as.character(Seurat::Idents(obj))
sce_out$seurat_clusters <- as.character(Seurat::Idents(obj))
zellkonverter::writeH5AD(sce_out, out)

metrics <- list(
  engine = annotation_method,
  phase = "downstream",
  resolution = res,
  n_clusters = length(unique(Seurat::Idents(obj))),
  n_cells = ncol(obj),
  scagent_annotation_method = annotation_method,
  azimuth_reference = if (is.null(azimuth_ref)) "none" else azimuth_ref,
  reference_only = TRUE
)
write(jsonlite::toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), file = "r_metrics.json")
cat("SCAGENT_R_OK annotate method=", annotation_method, " ref_only=TRUE\n", sep = "")
