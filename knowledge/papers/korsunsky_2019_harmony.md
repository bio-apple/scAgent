# Korsunsky et al. 2019 — Harmony: fast, sensitive and accurate integration

Nature Methods. 在 PCA 空间做迭代软聚类 + 簇内线性校正。

## 适用

- 输入是 PCA embedding，不是原始 counts。校正结果写入 `obsm['X_pca_harmony']`。
- 不改表达矩阵。下游邻居图 / UMAP / Leiden 用校正后的 embedding。
- 速度快，适合到百万细胞；可同时校正多个协变量（batch、donor、platform）。

## 限制

- 线性校正假设。强非线性批次或完全不共享细胞类型的数据集可能 overcorrect 或失败。
- 不能替代对实验设计的理解：处理效应如果与 batch 完全共线，Harmony 会把它当批次去掉。

## 与 scVI 的分工

- 要快、要可重复、批次结构简单 → Harmony。
- 要生成模型、缺失基因、多模态、不确定性 DE → scVI/scANVI（Lopez 2018; Gayoso 2022 scvi-tools）。
