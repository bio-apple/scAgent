# Squair et al. 2021 — Confronting false discoveries in single-cell differential expression

Nature Communications. 细胞级 DE 把细胞当独立样本，造成严重假阳性。

## 核心规则

- 生物学重复才是样本。同一供体的 5000 个 T 细胞不是 5000 个独立观测。
- 组间结论（疾病 vs 对照、处理 vs 对照）应 **pseudobulk**：按 sample × cell type 加和 counts，再用 DESeq2 / edgeR / limma-voom。
- Wilcoxon / t-test 在 cell-level 只适合探索 marker、看簇内富集，不能当论文主结论的 p 值。
- Crowell et al. 2020 (muscat) 给出了多样本、多细胞类型的 pseudobulk 框架。

## 智能体审查清单

- 代码是否按 `sample` 聚合？
- 是否报告了重复数 N？
- 是否做了多重校正（FDR / BH）？
- 单样本研究必须在报告里声明“探索性，无法做组间推断”。
