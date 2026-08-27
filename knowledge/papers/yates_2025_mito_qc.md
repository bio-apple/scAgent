# Yates et al. 2025 — Reconsidering mitochondrial QC cutoffs

Genome Biology. 质疑把高线粒体比例细胞一律当作濒死细胞过滤掉。

## 要点

- 线粒体比例升高可能是技术问题（破膜），也可能是真实生物学（代谢活跃、心肌、部分肿瘤、应激）。
- 固定 10% 或 5% cutoff 会在肿瘤和代谢研究中系统丢掉细胞。
- 正确做法：看 `pct_counts_mt` 与 `n_counts` 的 scatter；用 MAD 找离群；结合组织先验；在报告中写明移除数量与理由。

## 与 OSCA / emptyDrops 的关系

- Cell calling 用 barcode rank / emptyDrops，不要只靠 n_genes>200。
- QC 可视化三件套：**Violin、Scatter、MAD 判断**（本项目 qc_expert 硬性输出）。
