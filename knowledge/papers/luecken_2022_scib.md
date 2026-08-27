# Luecken et al. 2022 — scIB: Benchmarking atlas-level data integration

Nature Methods. 68 种整合方法在图谱尺度上的基准。

## 关键结论

- 没有万能整合器。评估必须同时看 **batch removal** 与 **bio-conservation**。
- HVG 选择通常提升整合效果；scaling 有时会偏向抹掉生物学差异。
- Harmony、scVI、Scanorama 等在不同数据形态上各有胜负。方法选择看任务：保留生物变异 vs 消除批次。
- UMAP 上“混匀”不是整合成功的证据。必须用 kBET、iLISI、cLISI、NMI/ARI 对标签、或 marker 可分性来检查 overcorrection。

## 对 scAgent 的规则

- 多样本时先问：批次是技术噪声还是生物学（处理、疾病、供体）？
- 默认 Python 路径：先 Harmony（快、可解释）；复杂非线性批次或需要不确定性时再用 scVI。
- 整合后必须审查：细胞类型是否被拉到一起但条件差异被抹平。

## 引用场景

用户问 “Harmony 还是 scVI”“要不要去批次”“整合效果怎么看”。
