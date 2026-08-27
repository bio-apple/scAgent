# Single-Cell RNA-seq Analysis Agent 单细胞生信分析智能体 

    scAgent/
    │
    ├── agents/                     # 智能体定义 (System Prompts, Tools)
    │   ├── planner/                # 任务拆解与全局调度
    │   ├── qc_expert/              # 动态质控策略制定
    │   ├── bio_coder/              # 代码生成 (支持 R/Seurat 与 Python/Scanpy)
    │   ├── annotation/             # 细胞注释与 Marker 校验
    │   ├── reviewer/               # 结果评估、异常检测与代码自纠错 (Debug)
    │   └── writer/                 # 整合图表与分析结论，生成报告
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
    

# skills 参考来源

SciAgent-Skills：https://github.com/jaechang-hits/SciAgent-Skills/tree/main/skills/genomics-bioinformatics/single-cell
