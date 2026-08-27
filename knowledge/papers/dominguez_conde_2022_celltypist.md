# Domínguez Conde et al. 2022 — Cross-tissue immune cell analysis (CellTypist)

Science. 逻辑回归参考注释；提供多组织免疫模型。

## 用法

- 输入：log1p 归一化到约 10^4 counts/cell 的 AnnData。
- 输出：per-cell 标签 + 可选 majority vote（按 Leiden 簇平滑）。
- 置信度 < 0.5 的细胞必须人工审查，可能是 doublet、过渡态或参考缺失类型。

## 分层证据（注释不允许单基因定论）

1. 无偏聚类（Leiden，resolution 由稳定性和 marker 可解释性校准）。
2. 参考映射（CellTypist / Azimuth / scANVI / popV）。
3. 至少两个独立 canonical marker + 阴性 marker 验证。
4. 与组织/疾病语境一致的命名（HCA 社区命名）。

## 失效

- 肿瘤微环境、发育阶段错配、跨物种直接套用成人免疫模型。
- popV（Luecken 2024）在需要 ensemble 不确定性时作为升级路径。
