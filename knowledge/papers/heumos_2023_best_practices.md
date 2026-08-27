# Heumos et al. 2023 — Best practices for single-cell analysis across modalities

Nature Reviews Genetics. 跨模态单细胞分析的共识工作流。

## 对分析智能体的约束

- 先看数据分布，再定 QC 阈值。不把 pctMT=10%、nHVG=2000 当铁律。
- 标准 scRNA-seq 路线：counts → QC（empty droplet / doublet / mito）→ 归一化 → HVG → 降维 → 邻域图 → 聚类 → 注释 → 差异表达。
- 组间差异表达应以生物学重复（样本）为单位；细胞不是独立样本。探索性 marker 可用 cell-level Wilcoxon，结论性比较用 pseudobulk + DESeq2/edgeR/limma。
- 整合不是默认步骤：构建参考图谱时需要；比较本就同质的样本时，先考虑 merge。
- 空间数据必须声明分辨率（spot / 亚细胞 / 单细胞）。spot 级结论不能写成单细胞事实。
- 报告必须包含软件版本、随机种子、过滤记录和参考数据集版本。

## 失效条件

- 极小样本或严重损坏的文库没有稳定分布，MAD 也可能失效。
- 参考偏倚会让注释的每一层都错（跨物种、跨发育阶段、肿瘤 vs 正常）。
