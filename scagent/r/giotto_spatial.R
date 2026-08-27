# Giotto spatial workflow stub (R)
# Usage: Rscript giotto_spatial.R <input_h5ad> <output_prefix>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript giotto_spatial.R input.h5ad output_prefix")
inp <- args[[1]]
prefix <- args[[2]]

if (!requireNamespace("Giotto", quietly = TRUE)) stop("Giotto not installed")
if (!requireNamespace("zellkonverter", quietly = TRUE)) stop("zellkonverter not installed")

sce <- zellkonverter::readH5AD(inp)
# Minimal placeholder: export coordinates if present for downstream Squidpy fallback
coords <- tryCatch({
  if ("spatial" %in% names(SingleCellExperiment::reducedDims(sce))) {
    as.data.frame(SingleCellExperiment::reducedDim(sce, "spatial"))
  } else {
    NULL
  }
}, error = function(e) NULL)

if (!is.null(coords)) write.csv(coords, paste0(prefix, "_spatial_coords.csv"), row.names = TRUE)

metrics <- list(engine = "giotto", phase = "spatial", note = "giotto stub; full pipeline may need Python Squidpy fallback")
write(jsonlite::toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), file = "r_metrics.json")
cat("SCAGENT_R_OK giotto_stub\n")
