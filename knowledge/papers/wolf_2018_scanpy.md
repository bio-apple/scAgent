# Wolf et al. 2018 — SCANPY: large-scale single-cell gene expression data analysis

Genome Biology. Python 可扩展 scRNA-seq 框架，AnnData 为中心对象。

## 与本仓库 skills 的对应

- 标准探索：`skills/scanpy-scrna-seq`
- 数据结构：`skills/anndata-data-structure`
- 聚类默认 Leiden（Traag 2019），不要无理由用 Louvain。
- 邻居图 → UMAP 只是可视化；聚类在邻域图上做，不在 UMAP 坐标上做。

## 平台提示

- 10x Genomics：`sc.read_10x_mtx` / `sc.read_10x_h5`，基因名 `var_names='gene_symbols'`。
- Parse Biosciences：注意不同 barcode 结构，planner 需标记 platform=parse，不要套 10x chemistry 假设。
- 多样本：先 concat 并保留 sample 列，再决定是否 Harmony/scVI。
