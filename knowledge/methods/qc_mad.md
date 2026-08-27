# QC 决策：数据驱动，而不是固定阈值

必须输出：

1. Violin：`n_genes_by_counts`、`total_counts`、`pct_counts_mt`（及 ribo/hb 如适用）
2. Scatter：`total_counts` vs `n_genes_by_counts`；`total_counts` vs `pct_counts_mt`
3. MAD：`median ± n * MAD`，n 通常 3–5，按组织调整

组织先验见 `config.yaml` 的 `qc_profiles`。心肌、肿瘤、代谢研究放宽线粒体过滤（Yates 2025）。
每步记录移除细胞数。空液滴用 barcode rank，不只 min_genes。
