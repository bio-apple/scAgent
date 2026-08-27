# Optional MAST cluster markers. Not a sample-level condition DE.
# CLI: Rscript mast.R counts.csv meta.csv out.csv
# counts.csv: genes x cells; meta.csv: cell,group  (cell = column names)

args <- commandArgs(trailingOnly = TRUE)
counts_path <- if (length(args) >= 1) args[[1]] else "mast_counts.csv"
meta_path <- if (length(args) >= 2) args[[2]] else "mast_meta.csv"
out_path <- if (length(args) >= 3) args[[3]] else "mast_de.csv"

.status <- function(status, detail = "") {
  txt <- sprintf('{"status":"%s","detail":"%s"}\n', status, gsub('"', "'", detail))
  writeLines(txt, out_path)
}

if (!requireNamespace("MAST", quietly = TRUE)) {
  .status("skipped", "MAST not installed")
  quit(save = "no", status = 0)
}

counts <- as.matrix(utils::read.csv(counts_path, row.names = 1, check.names = FALSE))
meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
if (!("cell" %in% names(meta)) || !("group" %in% names(meta))) {
  .status("skipped", "meta needs cell,group")
  quit(save = "no", status = 0)
}
cells <- as.character(meta$cell)
if (!all(cells %in% colnames(counts))) {
  .status("skipped", "meta$cell must match count columns")
  quit(save = "no", status = 0)
}
counts <- counts[, cells, drop = FALSE]
storage.mode(counts) <- "numeric"
mx <- max(counts, na.rm = TRUE)
if (is.finite(mx) && mx > 20) {
  counts <- log2(counts + 1)
}
sca <- MAST::FromMatrix(
  counts,
  cData = data.frame(wellKey = cells, group = factor(meta$group), row.names = cells),
  fData = data.frame(primerid = rownames(counts), row.names = rownames(counts))
)
levels <- sort(unique(as.character(meta$group)))
if (length(levels) != 2) {
  .status("skipped", "MAST helper needs exactly 2 groups")
  quit(save = "no", status = 0)
}
sca$group <- stats::relevel(factor(sca$group), ref = levels[[1]])
fit <- tryCatch(MAST::zlm(~group, sca), error = function(e) e)
if (inherits(fit, "error")) {
  .status("skipped", conditionMessage(fit))
  quit(save = "no", status = 0)
}
summ <- MAST::summary(fit, doLRT = paste0("group", levels[[2]]))
tab <- summ$datatable
hurdle <- tab[tab$component == "H" & tab$contrast == paste0("group", levels[[2]]), ]
logfc <- tab[tab$component == "logFC" & tab$contrast == paste0("group", levels[[2]]), c("primerid", "coef")]
names(logfc) <- c("gene", "logFC")
out <- merge(data.frame(gene = hurdle$primerid, pval = hurdle$`Pr(>Chisq)`), logfc, by = "gene", all.x = TRUE)
out$fdr <- stats::p.adjust(out$pval, method = "BH")
utils::write.csv(out[, c("gene", "logFC", "pval", "fdr")], out_path, row.names = FALSE)
cat("mast\n")
