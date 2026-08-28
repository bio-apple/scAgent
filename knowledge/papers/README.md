# knowledge/papers

把文献 PDF 或笔记 Markdown 放到本目录后执行：

```bash
python -m scagent ingest
```

智能体默认做融合检索（本目录 + `best_practices` / `methods` / `sops` / `upstream` / 结构化 `cell_ontology`·`marker_db`·`pathway` 等），作为 QC、整合、注释和统计审查的证据来源。
不要把受版权保护的全文提交到公开仓库；本地 PDF 仅用于私人 RAG 索引。
