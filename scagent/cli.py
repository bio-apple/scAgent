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
    hits = retrieve(args.query, collection=args.collection, top_k=args.top_k)
    print(format_hits(hits))
    return 0


def cmd_skills(_args: argparse.Namespace) -> int:
    for s in list_skills():
        print(f"{s.name}\n  {s.description}\n")
    return 0


def cmd_show_skill(args: argparse.Namespace) -> int:
    print(load_skill_text(args.name, include_references=args.refs))
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    from workflows.scRNA_langgraph import run_analysis

    cfg = load_config()
    ingest()
    _bm25_bundle.cache_clear()
    state = run_analysis(
        args.query,
        data_path=args.data,
        tissue=args.tissue,
        language=args.language or cfg["analysis"]["language"],
        execute_code=args.execute,
    )
    print(state.get("plan", {}).get("narrative", ""))
    print("\n--- skills ---")
    print(", ".join(state.get("skills_used") or []))
    print("\n--- review ---")
    print(json.dumps(state.get("review"), ensure_ascii=False, indent=2))
    print("\n报告: outputs/report.md")
    print("代码: workspace/analysis.py")
    return 0 if (state.get("review") or {}).get("passed", True) else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="scagent", description="单细胞生信分析智能体")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("ingest", help="索引 knowledge/papers（及 methods、markers）")
    s.set_defaults(func=cmd_ingest)

    s = sub.add_parser("retrieve", help="检索 RAG")
    s.add_argument("query")
    s.add_argument("--collection", default="papers")
    s.add_argument("--top-k", type=int, default=None)
    s.set_defaults(func=cmd_retrieve)

    s = sub.add_parser("skills", help="列出当前仓库已有 skills")
    s.set_defaults(func=cmd_skills)

    s = sub.add_parser("skill", help="打印某个 skill")
    s.add_argument("name")
    s.add_argument("--refs", action="store_true")
    s.set_defaults(func=cmd_show_skill)

    s = sub.add_parser("run", help="跑完整分析图（planner → QC → coder → review → 注释 → 报告）")
    s.add_argument("query", help="分析任务，例如：对 PBMC 做标准注释")
    s.add_argument("--data", default="", help="h5ad / 10x 目录")
    s.add_argument("--tissue", default="default")
    s.add_argument("--language", choices=["python", "r"], default=None)
    s.add_argument("--execute", action="store_true", help="在 workspace 中真正跑生成的脚本")
    s.set_defaults(func=cmd_run)
    return p


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
