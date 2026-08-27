# Pseudobulk DE: edgeR quasi-likelihood or DESeq2 on a gene x sample count matrix.
# CLI: Rscript deg.R counts.csv meta.csv out.csv [auto|edger|deseq2]
# counts.csv: row.names = gene, columns = sample ids (check.names=FALSE)
# meta.csv: sample, condition  (sample order must match columns)

run_deg <- function(counts_path, meta_path, out_path, engine = "auto") {
  counts <- as.matrix(utils::read.csv(counts_path, row.names = 1, check.names = FALSE))
  meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
  if (!("sample" %in% names(meta)) || !("condition" %in% names(meta))) {
    stop("meta.csv needs sample, condition")
  }
  samples <- as.character(meta$sample)
  if (!all(samples %in% colnames(counts))) {
    stop("meta$sample must match count columns")
  }
  counts <- counts[, samples, drop = FALSE]
  storage.mode(counts) <- "numeric"
  engine <- tolower(as.character(engine[[1]]))
  if (engine %in% c("auto", "edger", "edgeR")) {
    if (requireNamespace("edgeR", quietly = TRUE)) {
      res <- .edger_ql(counts, meta)
      utils::write.csv(res, out_path, row.names = FALSE)
      return("edger")
    }
    if (engine %in% c("edger", "edgeR")) {
      stop("edgeR is not installed")
    }
  }
  if (engine %in% c("auto", "deseq2")) {
    if (requireNamespace("DESeq2", quietly = TRUE)) {
      res <- .deseq2(counts, meta)
      utils::write.csv(res, out_path, row.names = FALSE)
      return("deseq2")
    }
    if (engine == "deseq2") {
      stop("DESeq2 is not installed")
    }
  }
  stop("no R DE backend (install edgeR and/or DESeq2)")
}

.edger_ql <- function(counts, meta) {
  group <- stats::relevel(factor(meta$condition), ref = sort(unique(as.character(meta$condition)))[[1]])
  y <- edgeR::DGEList(counts = counts, group = group)
  keep <- edgeR::filterByExpr(y, group = group)
  if (!any(keep)) {
    return(data.frame(gene = character(), logFC = numeric(), pval = numeric(), fdr = numeric()))
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  design <- stats::model.matrix(~group)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)
  tt <- as.data.frame(edgeR::topTags(qlf, n = Inf))
  data.frame(
    gene = rownames(tt),
    logFC = as.numeric(tt$logFC),
    pval = as.numeric(tt$PValue),
    fdr = as.numeric(tt$FDR),
    stringsAsFactors = FALSE
  )
}

.deseq2 <- function(counts, meta) {
  counts <- round(pmax(counts, 0))
  storage.mode(counts) <- "integer"
  colData <- data.frame(
    condition = stats::relevel(factor(meta$condition), ref = sort(unique(as.character(meta$condition)))[[1]]),
    row.names = as.character(meta$sample)
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~condition)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds)
  data.frame(
    gene = rownames(res),
    logFC = as.numeric(res$log2FoldChange),
    pval = as.numeric(res$pvalue),
    fdr = as.numeric(res$padj),
    stringsAsFactors = FALSE
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("usage: Rscript deg.R counts.csv meta.csv out.csv [auto|edger|deseq2]")
  }
  engine <- if (length(args) >= 4) args[[4]] else "auto"
  cat(run_deg(args[[1]], args[[2]], args[[3]], engine), sep = "\n")
}
