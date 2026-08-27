# Germain et al. 2020 / Scrublet — doublet detection before annotation

Doublet 会制造假细胞类型（两种 marker 的杂交）。注释前必须处理。

## 实践

- Scrublet / scDblFinder / SOLO 在 QC 阶段运行，分数写入 `obs`。
- 高 n_genes + 混合 marker 的簇优先怀疑 doublet。
- 不要把 doublet 簇注释成新细胞类型。
- 10x 高装载量文库 doublet 率更高，阈值随预期 doublet rate 调整。
