# Lopez et al. 2018 — Deep generative modeling for single-cell transcriptomics (scVI)

Nature Methods. 用 VAE 建模 UMI counts（ZINB/NB），潜空间分离生物学与批次。

## 实践要点

- 输入必须是 **raw counts**，不要 log-normalized 矩阵。counts 放 `layers['counts']`。
- `setup_anndata` → 实例化 `SCVI` → `train` → `get_latent_representation`。
- scANVI 用于半监督标签迁移；totalVI 用于 CITE-seq；DestVI 用于空间去卷积。
- 2025 年后的基础模型评测显示：简单基线在不少任务上不输 scFM。新模型必须在用户数据上小规模验证。

## 何时不要用

- <1 万细胞、批次简单、需要分钟级结果：Harmony 足够。
- 没有 GPU 且细胞数很大时训练成本高。
