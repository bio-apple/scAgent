# CellChat ligand-receptor analysis (R) — requires annotated Seurat
# Usage: Rscript cellchat.R <input_h5ad> <output_prefix> [groupby]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript cellchat.R input.h5ad output_prefix [groupby]")
inp <- args[[1]]
prefix <- args[[2]]
groupby <- if (length(args) >= 3) args[[3]] else "scagent_annotation"

suppressPackageStartupMessages({
  if (!requireNamespace("CellChat", quietly = TRUE)) stop("CellChat not installed")
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed")
  if (!requireNamespace("zellkonverter", quietly = TRUE)) stop("zellkonverter not installed")
})

sce <- zellkonverter::readH5AD(inp)
obj <- Seurat::CreateSeuratObject(counts = SummarizedExperiment::assay(sce, 1))
for (col in colnames(sce@colData)) obj[[col]] <- sce@colData[[col]]
if (!(groupby %in% colnames(obj@meta.data))) groupby <- "seurat_clusters"
Idents(obj) <- groupby

data.input <- Seurat::GetAssayData(obj, slot = "data")
meta <- obj@meta.data
cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = groupby)
CellChatDB <- CellChat::CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- CellChat::subsetData(cellchat)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::computeCommunProb(cellchat)
cellchat <- CellChat::aggregateNet(cellchat)
saveRDS(cellchat, paste0(prefix, "_cellchat.rds"))

metrics <- list(engine = "cellchat", phase = "cellchat", groupby = groupby)
write(jsonlite::toJSON(metrics, auto_unbox = TRUE, pretty = TRUE), file = "r_metrics.json")
cat("SCAGENT_R_OK cellchat\n")
