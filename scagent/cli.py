from __future__ import annotations

import argparse
import json
import sys

from rag.ingest import ingest
from rag.retriever import _bm25_bundle, format_hits, retrieve
from scagent.config import load_config
from scagent.skills_loader import list_skills, load_skill_text


def cmd_ingest(_args: argparse.Namespace) -> int:
    path = ingest()
    _bm25_bundle.cache_clear()
    print(f"indexed → {path}")
    return 0


def cmd_retrieve(args: argparse.Namespace) -> int:
    cols = args.collections.split(",") if args.collections else None
    hits = retrieve(args.query, collection=args.collection, top_k=args.top_k, collections=cols)
    print(format_hits(hits))
    return 0


def cmd_skills(_args: argparse.Namespace) -> int:
    for s in list_skills():
        print(f"{s.name}\n  {s.description}\n")
    return 0


def cmd_show_skill(args: argparse.Namespace) -> int:
    print(load_skill_text(args.name, include_references=args.refs))
    return 0


def cmd_demo(args: argparse.Namespace) -> int:
    from scagent.demo import DEFAULT_PATH, write_tiny_h5ad

    path = write_tiny_h5ad(args.path or DEFAULT_PATH)
    print(path)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    from workflows.scRNA_langgraph import run_analysis

    cfg = load_config()
    ingest()
    _bm25_bundle.cache_clear()
    mode = "full"
    if args.qc_only:
        mode = "qc_only"
    elif args.annotate_only:
        mode = "annotate_only"
    execute = False if args.dry_run else args.execute
    state = run_analysis(
        args.query,
        data_path=args.data,
        tissue=args.tissue,
        language=args.language or cfg["analysis"]["language"],
        execute_code=execute,
        mode=mode,
        interrupt_after_qc=args.interrupt,
        auto_confirm=args.yes or not args.interrupt,
        resolution=args.resolution,
        batch_key=args.batch_key,
        markers_path=args.markers,
        report_lang=args.report_lang,
        integrator=args.integrator,
        imputation=args.impute,
        qc_method=args.qc_method,
    )
    print(state.get("plan", {}).get("narrative", ""))
    print("\n--- status ---", state.get("status"))
    print("\n--- skills ---")
    print(", ".join(state.get("skills_used") or []))
    print("\n--- review QC ---")
    print(json.dumps(state.get("review_qc"), ensure_ascii=False, indent=2))
    print("\n--- review downstream ---")
    print(json.dumps(state.get("review_downstream"), ensure_ascii=False, indent=2))
    print("\n报告: outputs/report.md")
    print("脚本: workspace/qc_preprocess.py  workspace/cluster_annotate.py  workspace/reproducible_script.py")
    if state.get("status") == "awaiting_qc_confirmation":
        print("QC 已暂停。检查 workspace/ 后用 --annotate-only --yes 继续。")
    qc_ok = (state.get("review_qc") or {}).get("passed", True)
    down_ok = (state.get("review_downstream") or {}).get("passed", True)
    if state.get("r_degraded"):
        return 2
    return 0 if qc_ok and down_ok else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="scagent", description="单细胞生信分析智能体")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("ingest", help="索引 knowledge/papers、methods、markers 与 best_practices/reference")
    s.set_defaults(func=cmd_ingest)

    s = sub.add_parser("retrieve", help="检索 RAG")
    s.add_argument("query")
    s.add_argument("--collection", default="papers")
    s.add_argument("--collections", default=None, help="逗号分隔，如 papers,markers,best_practices")
    s.add_argument("--top-k", type=int, default=None)
    s.set_defaults(func=cmd_retrieve)

    s = sub.add_parser("skills", help="列出当前仓库已有 skills")
    s.set_defaults(func=cmd_skills)

    s = sub.add_parser("skill", help="打印某个 skill")
    s.add_argument("name")
    s.add_argument("--refs", action="store_true")
    s.set_defaults(func=cmd_show_skill)

    s = sub.add_parser("demo", help="写出 100 细胞稀疏 demo h5ad")
    s.add_argument("--path", default=None)
    s.set_defaults(func=cmd_demo)

    s = sub.add_parser("run", help="两阶段图：QC 执行审查 → 聚类注释审查 → 报告")
    s.add_argument("query", help="分析任务，例如：对 PBMC 做标准注释")
    s.add_argument("--data", default="", help="h5ad / 10x 目录 / Seurat .rds")
    s.add_argument("--tissue", default="default")
    s.add_argument("--language", choices=["python", "r"], default=None)
    s.add_argument("--execute", action="store_true", help="在 workspace 中真正跑生成的脚本")
    s.add_argument("--dry-run", action="store_true", help="只生成脚本与报告，不执行")
    s.add_argument("--qc-only", action="store_true")
    s.add_argument("--annotate-only", action="store_true", help="跳过 QC，需已有 workspace/adata_qc.h5ad")
    s.add_argument("--interrupt", action="store_true", help="QC 审查通过后暂停，供人工确认阈值")
    s.add_argument("--yes", action="store_true", help="跳过人工确认")
    s.add_argument("--resolution", type=float, default=None)
    s.add_argument("--batch-key", default=None)
    s.add_argument("--integrator", choices=["auto", "none", "harmony", "scvi", "cca"], default=None, help="批次校正模块，默认读 config.modules.batch")
    s.add_argument("--impute", choices=["none", "magic", "alra"], default=None, help="Dropout 插补，默认读 config.modules.imputation")
    s.add_argument("--qc-method", choices=["mad", "percentile", "hybrid"], default=None)
    s.add_argument("--markers", default=None, help="自定义 marker CSV/JSON（可含 lineage 列）")
    s.add_argument("--report-lang", choices=["zh", "en", "both"], default="zh")
    s.set_defaults(func=cmd_run)
    return p


def main(argv: list[str] | None = None) -> None:
    load_config()
    parser = build_parser()
    args = parser.parse_args(argv)
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
