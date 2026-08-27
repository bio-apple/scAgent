# Single-Cell RNA-seq Analysis Agent 单细胞生信分析智能体 

    scAgent/
    │
    ├── agents/                     # 智能体定义 (System Prompts, Tools)
    │   ├── planner/                # 调度核心:读取 metadata、判断物种、平台（10x/Parse）、是否多样本，然后决定分析路线
    │   ├── qc_expert/              # 动态质控策略制定:不直接用固定阈值，而是根据组织自动推荐 QC，例如肿瘤、脑组织、PBMC 不同策略。输出必须包含 VlnPlot、Scatter、MAD 判断。
    │   ├── bio_coder/              # 代码生成 主语言 R（Seurat)，Python:Scanpy、Squidpy（备用）
    │   ├── annotation/             # 细胞注释与 Marker 校验,综合 Marker + CellMarker2.0 + PanglaoDB + Azimuth，不允许只看一个基因决定细胞类型。
    │   ├── reviewer/               # 结果评估、异常检测与代码自纠错 (Debug),审核每一步是否符合统计规范，例如 DEG 是否做多重校正、批次是否真的消除、UMAP 是否过聚类。
    │   └── writer/                 # 整合图表与分析结论，生成报告,自动写 Result，不解释不存在的现象，每张图都有 Figure legend。
    │
    ├── skills/                     # 标准操作规程 (SOP / Standard Prompting)
    │   ├── R/
    │   │   ├── seurat_qc.md
    │   │   ├── harmony.md
    │   │   └── cellchat.md
    │   └── python/
    │       ├── scanpy_qc.md
    │       └── scrublet.md
    │
    ├── knowledge/                  # RAG 知识库
    │   ├── methods/                # 算法原理与适用场景 (如 Scran vs LogNormalize)
    │   ├── markers/                # CellMarker, PanglaoDB 数据库
    │   └── papers/                 # 最新文献洞察
    │
    ├── workflows/                  # 基于 LangGraph 的状态图定义
    │   ├── state.py                # 定义 AgentState (传递全局分析数据与代码状态)
    │   └── scRNA_langgraph.py      # 主流程编排
    │
    ├── sandbox/                    # 代码执行沙盒 (Jupyter Kernel / Docker)
    │   └── executor.py
    │
    ├── workspace/                  # 运行时动态生成的临时代码与图像
    ├── outputs/                    # 终版产物 (RDS/H5AD/HTML Report)
    │
    ├── prompts/                    # 公用 Prompt 模版
    ├── report_templates/           # RMarkdown / Quarto / HTML 报告模板
    ├── tests/                      # 单元测试 (测试 Agent 是否能正确纠错)
    ├── requirements.txt            # Python 依赖
    └── config.yaml                 # 模型 API Key、内存限制、路径配置

# 技术栈

|层 | 推荐                 |
|-|--------------------|
|Agent | LangGraph          |




# skills 参考来源

SciAgent-Skills：https://github.com/jaechang-hits/SciAgent-Skills/tree/main/skills/genomics-bioinformatics/single-cell
