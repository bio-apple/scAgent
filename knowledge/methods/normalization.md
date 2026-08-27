# 归一化：LogNormalize vs scran vs SCTransform

- **LogNormalize**（Seurat / Scanpy `normalize_total` + `log1p`）：简单、可复现，大多数探索性分析够用。
- **scran deconvolution**：用池化细胞估计 size factor，对组成差异更稳健，偏统计传统（Lun / Bioconductor）。
- **SCTransform**：把技术噪声作为正则化负二项回归；叙事不同，与 log 结果高度相关。

智能体默认 Python 路径用 `normalize_total(target_sum=1e4)` + `log1p`，与 CellTypist 输入约定一致。若用户明确要 scran/SCT，再切换。
