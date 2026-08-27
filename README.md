# Single-Cell RNA-seq Analysis Agent

基于 LangGraph 的单细胞生信分析智能体：调度 → 动态 QC → 代码生成 → 审查纠错 → 注释 → 报告。  
RAG 默认检索 `knowledge/papers`。仓库里现有的 SciAgent-style skills **保持原样**，作为可执行 SOP。

## 目录

```
scAgent/
├── agents/                     # 智能体
│   ├── planner.py              # 读 metadata，判断物种 / 平台（10x/Parse）/ 多样本，决定路线
│   ├── qc_expert.py            # 组织感知 QC；必须含 Violin、Scatter、MAD
│   ├── bio_coder.py            # 主语言 Python/Scanpy（对齐现有 skills）；R/Seurat 为显式备选
│   ├── annotation.py           # Marker + 参考映射，禁止单基因定论
│   ├── reviewer.py             # 统计规范、过聚类、假整合、DEG 多重校正
│   └── writer.py               # 报告；不解释不存在的现象；图注
├── skills/                     # 已有 SOP（不要拆到 R/python 子目录）
│   ├── scanpy-scrna-seq/
│   ├── anndata-data-structure/
│   ├── harmony-batch-correction/
│   ├── scvi-tools-single-cell/
│   ├── celltypist-cell-annotation/
│   ├── popv-cell-annotation/
│   ├── cellxgene-census/
│   └── single-cell-annotation-guide/
├── knowledge/                  # RAG 语料
│   ├── papers/                 # 默认检索集合（可再放入 PDF）
│   ├── methods/
│   └── markers/
├── workflows/
│   ├── state.py
│   └── scRNA_langgraph.py
├── sandbox/executor.py
├── prompts/
├── report_templates/
├── tests/
├── config.yaml
└── requirements.txt
```

流程：

`inspect → planner → qc_expert → bio_coder → execute → reviewer ⇄ bio_coder → annotation → writer`

## 安装

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

可选：分析栈 `pip install -e ".[analysis]"`。  
LLM 可选。未配置 `OPENAI_API_KEY` 时走确定性规划 + Scanpy 模板，图仍可跑通。

```bash
cp .env.example .env   # 填写 OPENAI_API_KEY；兼容 OpenAI 接口可设 OPENAI_BASE_URL
```

## 使用

```bash
# 索引 knowledge/papers
python -m scagent ingest

# 检索
python -m scagent retrieve "Harmony versus scVI"

# 已有 skills
python -m scagent skills

# 规划 + 生成脚本 + 报告（默认不跑分析）
python -m scagent run "对 PBMC 做标准分析并注释" --data /path/to/data.h5ad --tissue pbmc

# 在 workspace/ 真正执行生成的代码
python -m scagent run "..." --data /path/to/data.h5ad --execute
```

产物：`workspace/analysis.py`、`outputs/report.md`。把 PDF 放进 `knowledge/papers/` 后重新 `ingest` 即可进入 RAG。

## 设计选择

- **Skills**：保留当前 `skills/*/SKILL.md`，不按 README 草稿改成 `skills/R` 与 `skills/python`。
- **语言**：现有 skills 都是 Scanpy 生态，因此 bio_coder 默认 Python；用户传 `--language r` 时会警告缺少对等 SOP。
- **RAG**：BM25 检索 `knowledge/papers`（同时入库 methods/markers）。不依赖外部 embedding 服务。
- **QC**：组织 profile 在 `config.yaml`；硬性三件套 Violin / Scatter / MAD。
- **审查**：缺 QC 三件套或多样本无整合且无理由 → 打回 bio_coder（最多 2 次）。

## 测试

```bash
pytest -q
```

Skills 参考来源：[SciAgent-Skills single-cell](https://github.com/jaechang-hits/SciAgent-Skills/tree/main/skills/genomics-bioinformatics/single-cell)
