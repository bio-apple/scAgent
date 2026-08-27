# Traag, Waltman & van Eck 2019 — From Louvain to Leiden

Scientific Reports. Leiden 修复 Louvain 的 disconnected community 问题。

## 对聚类的约束

- 默认 Leiden。resolution 不是 0.8 教条，要用簇稳定性、marker 可分性、生物学粒度来校准。
- 过聚类：同一细胞类型碎成许多簇 → reviewer 应警告 UMAP/Leiden 过细。
- 欠聚类：明显不同 lineage 被合并 → 提高 resolution 或分层聚类。
- 发育连续谱上，离散簇是近似；需要轨迹时改用 PAGA / 拟时序，并声明假设。
