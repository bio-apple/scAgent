# knowledge/ — RAG + 结构化知识库

所有可检索文档都在本目录。Planner / QC / 注释 / 解读走 **融合检索**；细胞类型、通路、疾病签名走 **结构化 JSON 记录**（`scagent.kb.lookup_structured`），不是 prompt 里的散文。

| 子目录 | collection | 来源 | 用途 |
|--------|------------|------|------|
| `cell_ontology/` | `cell_ontology` | Cell Ontology (CL) 子集 | 细胞类型 ID / 同义词 |
| `marker_db/` | `marker_db` | CellMarker 2.0 / PanglaoDB 子集 | 双阳性 + 阴性 marker |
| `pathway/` | `pathway` | MSigDB Hallmark / GO 子集 | 离线 ORA 基因集 |
| `disease_signature/` | `disease_signature` | 文献签名（≥2 marker + GO + DOI） | 细胞状态证据链 |
| `tissue_reference/` | `tissue_reference` | Human Cell Atlas 器官图谱子集 | 组织预期细胞类型 |
| `best_practices/` | `best_practices` | Heumos 2023 步骤 SOP | 分析路线 |
| `papers/` | `papers` | 文献笔记 / 本地 PDF | 证据 |
| `methods/` | `methods` | scAgent 方法卡 | QC / 整合决策 |
| `sops/` | `sops` | 实验室 SOP | `add-doc` |
| `upstream/` | `upstream` | theislab 全书 | `update-kb`（gitignore） |
| `evidence/` | — | 兼容旧路径 | 现指向 `disease_signature/` |
| `markers/` | — | 兼容旧路径 | 现指向 `marker_db/` |

以上数据库均为 **离线精选子集**，不是全库镜像。

```bash
python -m scagent ingest
python -m scagent retrieve "CL:0000084 T cell"
python -m scagent retrieve "HALLMARK_HYPOXIA" --collections pathway,papers
```
