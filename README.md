# Single-Cell RNA-seq Analysis Agent

> **Languages:** 中文（本文） | [English](README.en.md)

LangGraph 单细胞生信智能体：**组织感知 QC、执行审查、可审计脚本**。内置 99 个 skills + [`knowledge/`](README.en.md#knowledge-base) 融合 RAG / 结构化 KB。

## 快速开始

```bash
pip install -r requirements.txt && pip install -e ".[dev]"
python -m scagent init
python -m scagent run "PBMC 标准分析" --data tests/data/tiny_100cells.h5ad --tissue pbmc --dry-run
```

真跑：`pip install -r requirements-analysis.txt`，加 `--execute`。API key 用环境变量 `OPENAI_API_KEY`（可选，无 key 为模板模式）。

```bash
python -m scagent ingest                    # 首次或更新知识库后
python -m scagent run "…" --data data.h5ad --tissue pbmc --execute
python -m scagent run --help                  # 全部 CLI 参数
```

**完整文档**（架构、CLI、知识库、Tool Router、策略与排错）：[README.en.md](README.en.md)

## 常见问题

| 现象 | 处理 |
|------|------|
| 报告「未执行」 | `--execute` + `requirements-analysis.txt` |
| `retrieve` 无结果 | `scagent ingest` |
| `--resume` 版本冲突 | `--force-resume` 或新开 run |

```bash
pytest -q    # CI: pytest + flake8 + black
```
