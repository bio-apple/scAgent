from __future__ import annotations

import argparse
import json
import sys

from rag.ingest import ingest
from rag.retriever import clear_retrieve_cache, format_hits, retrieve
from scagent.config import load_config
from scagent.skills_loader import list_skills, load_skill_text


def cmd_ingest(_args: argparse.Namespace) -> int:
    path = ingest(force=True)
    clear_retrieve_cache()
    print(f"indexed → {path}")
    return 0


def cmd_update_kb(args: argparse.Namespace) -> int:
    from rag.kb import update_kb

    try:
        info = update_kb(url=args.url, branch=args.branch)
    except Exception as exc:
        print(f"update-kb 失败: {exc}")
        return 1
    commit = info.get("commit") or "unknown"
    print(f"sc-best-practices → {info['dest']} ({info['n_files']} files, {commit})")
    if info.get("index"):
        print(f"indexed → {info['index']}")
    if info.get("n_files") == 0:
        print("未找到 jupyter-book 章节（md/ipynb）。检查 --url / --branch。")
        return 1
    return 0


def cmd_add_doc(args: argparse.Namespace) -> int:
    from rag.kb import add_doc

    try:
        info = add_doc(args.path, name=args.name)
    except Exception as exc:
        print(f"add-doc 失败: {exc}")
        return 1
    print(f"added {info['n_files']} file(s) → {info['dest']}")
    for p in info.get("files") or []:
        print(f"  {p}")
    if info.get("index"):
        print(f"indexed → {info['index']}")
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


def cmd_snapshots(args: argparse.Namespace) -> int:
    from scagent.snapshot import list_snapshots
    from workflows.checkpointing import load_last_thread

    tid = args.thread_id or load_last_thread()
    if not tid:
        print("无 thread_id（先 run 或传 --thread-id）")
        return 1
    entries = list_snapshots(tid)
    if not entries:
        print(f"thread {tid}: 尚无 AnnData 快照")
        return 0
    print(f"thread {tid}")
    for e in entries:
        print(f"  {e.get('step')}\tkind={e.get('kind')}\tstored={e.get('stored_bytes')}\t{e.get('path')}")
    return 0


def cmd_branch(args: argparse.Namespace) -> int:
    from scagent.config import load_config, resolve_path
    from scagent.snapshot import checkout, fork_branch
    from workflows.checkpointing import load_last_thread

    src = args.from_thread or load_last_thread()
    if not src:
        print("需要 --from-thread 或已有 last_thread_id")
        return 1
    payload = fork_branch(src, args.as_name, from_step=args.step)
    print(f"branch {src} → {args.as_name} steps={[e.get('step') for e in payload.get('entries') or []]}")
    if args.checkout:
        ws = resolve_path(load_config(), "workspace")
        dest = checkout(args.as_name, args.step or "qc", ws)
        print(f"checkout → {dest}")
    return 0


def cmd_view(args: argparse.Namespace) -> int:
    from scagent.config import load_config, resolve_path
    from scagent.viewer import export_workspace_viewer, serve_viewer

    out = resolve_path(load_config(), "outputs")
    h5ad = args.data or None
    if args.serve:
        serve_viewer(out, port=args.port, h5ad=h5ad)
        return 0
    path = export_workspace_viewer(out, h5ad=h5ad)
    if path is None:
        print("没有 adata_processed.h5ad / adata_qc.h5ad。先 --execute 或指定 --data")
        return 1
    print(path)
    if args.open:
        import webbrowser

        webbrowser.open(path.resolve().as_uri())
    return 0


def cmd_ask(args: argparse.Namespace) -> int:
    from scagent.config import load_config, resolve_path
    from scagent.io import read_h5ad
    from scagent.viewer import find_workspace_h5ad, load_selection, summarize_selection

    sel = load_selection(args.selection)
    query = args.query or sel.get("query") or "分析我框选的这组细胞"
    h5ad = args.data or find_workspace_h5ad()
    if h5ad is None:
        print("没有 h5ad。用 --data 或先跑分析")
        return 1
    adata = read_h5ad(h5ad, backed=True)
    summary = summarize_selection(adata, sel["cell_ids"], query=query)
    out = resolve_path(load_config(), "outputs")
    out.mkdir(parents=True, exist_ok=True)
    (out / "selection.json").write_text(
        json.dumps({**sel, "query": query}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (out / "selection_report.md").write_text(summary["text"] + "\n", encoding="utf-8")
    print(summary["text"])
    print("报告: outputs/selection_report.md")
    return 0 if summary["n_matched"] else 1


def cmd_memory(_args: argparse.Namespace) -> int:
    from agents.memory import dump_memory_yaml, load_memory

    mem = load_memory()
    if not mem:
        print("尚无 memory.yaml（先跑 scagent run）")
        return 1
    print(dump_memory_yaml(mem).rstrip())
    return 0


def cmd_demo(args: argparse.Namespace) -> int:
    from scagent.demo import DEFAULT_PATH, write_tiny_h5ad

    path = write_tiny_h5ad(args.path or DEFAULT_PATH)
    print(path)
    return 0


def cmd_confirm(args: argparse.Namespace) -> int:
    from scagent.config import load_config, resolve_path
    from scagent.hitl import load_session, save_session
    from workflows.scRNA_langgraph import run_analysis

    session = load_session()
    if not (session.get("user_query") or "").strip():
        print("没有 HITL 会话。先: scagent run \"...\" --data ... --interrupt")
        return 1
    qc_choice = session.get("qc_choice")
    resolution_choice = session.get("resolution_choice")
    if args.kind == "mt":
        qc_choice = args.choice
    else:
        resolution_choice = args.choice
        if not qc_choice:
            qc_choice = "recommended"
    save_session({**session, "qc_choice": qc_choice, "resolution_choice": resolution_choice})
    mode = "full"
    ws = resolve_path(load_config(), "workspace")
    if args.kind == "resolution" and (ws / "adata_qc.h5ad").exists():
        mode = "annotate_only"
    execute = False if args.dry_run else bool(args.execute)
    state = run_analysis(
        session["user_query"],
        data_path=session.get("data_path") or "",
        tissue=session.get("tissue") or "default",
        language=session.get("language"),
        execute_code=execute,
        mode=mode,
        interrupt_after_qc=True,
        auto_confirm=False,
        qc_choice=qc_choice,
        resolution_choice=resolution_choice,
        markers_path=session.get("markers_path"),
        batch_key=session.get("batch_key"),
        report_lang=session.get("report_lang") or "zh",
        thread_id=args.thread_id or session.get("thread_id"),
    )
    print("\n--- status ---", state.get("status"), "thread_id=", state.get("thread_id"))
    if state.get("status") == "awaiting_mt_confirmation":
        print("继续: scagent confirm mt recommended|lenient|strict")
    elif state.get("status") == "awaiting_resolution_confirmation":
        print("继续: scagent confirm resolution recommended|coarse|fine")
    print("决策卡: outputs/decisions/mt.html  outputs/decisions/resolution.html")
    return 0


def cmd_init(args: argparse.Namespace) -> int:
    from scagent.init_wizard import (
        _as_bool,
        build_run_argv,
        collect_session,
        format_command,
        overlay_path,
        resource_overlay,
        session_path,
        summarize,
        write_overlay,
        write_session,
    )

    cfg = load_config()
    answers = {
        "data": args.data,
        "tissue": args.tissue,
        "task": args.task,
        "query": args.query,
        "language": args.language,
        "report_lang": args.report_lang,
        "memory_mb": args.memory_mb,
        "n_jobs": args.n_jobs,
        "timeout": args.timeout,
        "batch_key": args.batch_key,
        "condition_key": args.condition_key,
    }
    if args.execute:
        answers["execute"] = True
    elif args.dry_run:
        answers["execute"] = False
    if args.interrupt:
        answers["interrupt"] = True
    try:
        session = collect_session(answers=answers, use_defaults=bool(args.yes), cfg=cfg)
    except (EOFError, KeyboardInterrupt):
        print("已取消")
        return 1
    print()
    print(summarize(session))
    if not args.yes:
        try:
            ok = input("按上述配置写入？ [Y/n]: ")
        except (EOFError, KeyboardInterrupt):
            print("已取消")
            return 1
        if str(ok).strip().lower() in {"n", "no", "否"}:
            print("已取消")
            return 1
    sess_file = write_session(session, session_path(cfg))
    print(f"会话 → {sess_file}")
    overlay = resource_overlay(session, cfg)
    local = None
    if overlay and not args.no_write_config:
        local = write_overlay(overlay, overlay_path(cfg))
        if local:
            print(f"资源覆盖 → {local}（不改 config.yaml；密钥仍用环境变量）")
            load_config(reload=True)
    cmd = format_command(session)
    print("下次运行:")
    print(f"  {cmd}")
    start = bool(args.run)
    if start and not args.yes:
        try:
            start = _as_bool(input("现在开始分析？ [n]: "), False)
        except (EOFError, KeyboardInterrupt):
            start = False
    if not start:
        return 0
    return cmd_run(build_parser().parse_args(build_run_argv(session)))


def cmd_run(args: argparse.Namespace) -> int:
    from workflows.scRNA_langgraph import run_analysis
    from scagent.compat import ResumeIncompatibleError
    from scagent.viewer import load_selection

    cfg = load_config()
    ingest()
    clear_retrieve_cache()
    mode = "full"
    if args.qc_only:
        mode = "qc_only"
    elif args.annotate_only:
        mode = "annotate_only"
    execute = False if args.dry_run else args.execute
    selection = load_selection(args.selection) if getattr(args, "selection", None) else None
    query = args.query
    if selection and selection.get("query") and not query:
        query = selection["query"]
    try:
        state = run_analysis(
            query,
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
            remove_doublets=True if args.remove_doublets else None,
            doublet_methods=args.doublet_methods,
            ambient=args.ambient,
            condition_key=args.condition_key,
            deg_engine=getattr(args, "deg_engine", None),
            marker_method=getattr(args, "marker_method", None),
            deg_cross_validate=getattr(args, "deg_cross_validate", None),
            thread_id=args.thread_id,
            resume=args.resume,
            force_resume=bool(getattr(args, "force_resume", False)),
            selection=selection,
        )
    except ResumeIncompatibleError as exc:
        print(str(exc))
        print("版本不兼容时拒绝续跑。确认后可加 --force-resume。")
        return 2
    print(state.get("plan", {}).get("narrative", ""))
    print("\n--- status ---", state.get("status"), "thread_id=", state.get("thread_id"))
    print("\n--- skills ---")
    print(", ".join(state.get("skills_used") or []))
    from agents.reviewer import format_review_card

    print("\n--- Reviewer ---")
    print(format_review_card(state.get("review_publication"), args.report_lang).rstrip())
    print("\n--- review QC ---")
    print(json.dumps(state.get("review_qc"), ensure_ascii=False, indent=2))
    print("\n--- review downstream ---")
    print(json.dumps(state.get("review_downstream"), ensure_ascii=False, indent=2))
    from scagent.dual import render_dual_console

    print("\n" + render_dual_console(state).rstrip())
    print("\n报告: outputs/report.md")
    print("双重输出: outputs/dual.md")
    print("provenance: outputs/memory.yaml")
    print("notebook: outputs/analysis.ipynb （R 路径为 outputs/analysis.Rmd）")
    print("viewer: outputs/viewer.html （Plotly 框选细胞；scagent view --serve 可当场提问）")
    print("脚本: workspace/qc_preprocess.py  workspace/cluster_annotate.py  workspace/reproducible_script.py")
    if state.get("status") == "awaiting_mt_confirmation":
        print("线粒体阈值待确认。打开 outputs/decisions/mt.html")
        print("继续: scagent confirm mt recommended|lenient|strict")
    elif state.get("status") == "awaiting_resolution_confirmation":
        print("Leiden resolution 待确认。打开 outputs/decisions/resolution.html")
        print("继续: scagent confirm resolution recommended|coarse|fine")
    elif state.get("status") == "awaiting_qc_confirmation":
        print("QC 已暂停。检查 workspace/ 后用 scagent confirm 继续。")
    print("决策卡: outputs/decisions/mt.html  outputs/decisions/resolution.html")
    qc_ok = (state.get("review_qc") or {}).get("passed", True)
    down_ok = (state.get("review_downstream") or {}).get("passed", True)
    if state.get("r_degraded"):
        return 2
    return 0 if qc_ok and down_ok else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="scagent", description="单细胞生信分析智能体")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("init", help="交互式配置向导（数据路径、组织、任务、资源限制）")
    s.add_argument("--data", default=None, help="h5ad / loom / Cell Ranger outs / 10x mtx / Seurat .rds；逗号分隔或样本父目录")
    s.add_argument("--tissue", default=None)
    s.add_argument("--task", choices=["annotate", "qc", "deg", "trajectory", "custom"], default=None)
    s.add_argument("--query", default=None, help="任务描述，默认由任务类型生成")
    s.add_argument("--execute", action="store_true", help="真正执行生成的脚本")
    s.add_argument("--dry-run", action="store_true", help="只生成脚本（默认）")
    s.add_argument("--yes", action="store_true", help="其余项用默认值，不提问")
    s.add_argument("--run", action="store_true", help="写完配置后立刻 scagent run")
    s.add_argument("--language", choices=["python", "r"], default=None)
    s.add_argument("--report-lang", choices=["zh", "en", "both"], default=None)
    s.add_argument("--memory-mb", type=int, default=None, help="sandbox 内存上限 MB")
    s.add_argument("--n-jobs", type=int, default=None, help="Scanpy/joblib 并行度，-1=全部 CPU")
    s.add_argument("--timeout", type=int, default=None, help="单阶段超时秒")
    s.add_argument("--batch-key", default=None)
    s.add_argument("--condition-key", default=None)
    s.add_argument("--interrupt", action="store_true")
    s.add_argument("--no-write-config", action="store_true", help="不写 config.local.yaml")
    s.set_defaults(func=cmd_init)

    s = sub.add_parser("ingest", help="索引 knowledge/papers、methods、markers、best_practices 与 knowledge/sops")
    s.set_defaults(func=cmd_ingest)

    s = sub.add_parser("update-kb", help="一键拉取最新 sc-best-practices（theislab）并重建 RAG 索引")
    s.add_argument("--url", default=None, help="git URL，默认 https://github.com/theislab/single-cell-best-practices.git")
    s.add_argument("--branch", default="main", help="分支，默认 main（失败时 clone 会再试 master）")
    s.set_defaults(func=cmd_update_kb)

    s = sub.add_parser("add-doc", help="把实验室 SOP 纳入本地 RAG（复制到 knowledge/sops 后 ingest）")
    s.add_argument("path", help="文件或目录（md / txt / pdf / ipynb）")
    s.add_argument("--name", default=None, help="写入 sops 时的文件名或子目录名")
    s.set_defaults(func=cmd_add_doc)

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

    s = sub.add_parser("memory", help="打印分析 provenance（步骤与参数，不是聊天）")
    s.set_defaults(func=cmd_memory)

    s = sub.add_parser("snapshots", help="列出当前 thread 的惰性 AnnData 快照（硬链接/obs 增量，不把 X 拷进 RAM）")
    s.add_argument("--thread-id", default=None)
    s.set_defaults(func=cmd_snapshots)

    s = sub.add_parser("branch", help="从某步快照分叉参数实验：共享 h5ad inode，不复制矩阵")
    s.add_argument("--from-thread", default=None, help="源 thread_id，默认 last_thread_id")
    s.add_argument("--as", dest="as_name", required=True, help="新 branch / thread 名")
    s.add_argument("--step", default="qc", help="分叉点：qc 或 downstream")
    s.add_argument("--checkout", action="store_true", help="硬链接到 workspace/ 以便接着跑")
    s.set_defaults(func=cmd_branch)

    s = sub.add_parser("run", help="Planner → Executor → Reviewer → Publication Report")
    s.add_argument("query", help="分析任务，例如：对 PBMC 做标准注释")
    s.add_argument("--data", default="", help="h5ad / loom / Cell Ranger outs / 10x mtx / Seurat .rds；逗号分隔=多样本拼接")
    s.add_argument("--tissue", default="default")
    s.add_argument("--language", choices=["python", "r"], default=None)
    s.add_argument("--execute", action="store_true", help="在 workspace 中真正跑生成的脚本")
    s.add_argument("--dry-run", action="store_true", help="只生成脚本与报告，不执行")
    s.add_argument("--qc-only", action="store_true")
    s.add_argument("--annotate-only", action="store_true", help="跳过 QC，需已有 workspace/adata_qc.h5ad")
    s.add_argument("--interrupt", action="store_true", help="在线粒体阈值与 Leiden resolution 两处暂停，等人看直方图再确认")
    s.add_argument("--yes", action="store_true", help="跳过人工确认")
    s.add_argument("--resolution", type=float, default=None)
    s.add_argument("--batch-key", default=None)
    s.add_argument("--integrator", choices=["auto", "none", "harmony", "scvi", "cca", "scanorama", "bbknn"], default=None, help="批次校正。auto：inspect 检测到批次后 Harmony 或 scVI；cca/scanorama=Scanorama；bbknn 改邻居图")
    s.add_argument("--impute", choices=["none", "magic", "alra"], default=None, help="Dropout 插补，默认读 config.modules.imputation")
    s.add_argument("--ambient", choices=["auto", "none", "soupx", "decontx"], default=None, help="Ambient RNA 去除，默认读 config.modules.ambient")
    s.add_argument("--remove-doublets", action="store_true", help="按两法共识 predicted_doublet 过滤细胞")
    s.add_argument(
        "--doublet-methods",
        choices=["auto", "scrublet", "both"],
        default=None,
        help="auto：多样本/复杂组织 Scrublet+scDblFinder（无 R 则表达模拟）；both 强制交叉验证",
    )
    s.add_argument("--condition-key", default=None, help="组间比较的 obs 列名；出现时走 sample-level pseudobulk DE")
    s.add_argument(
        "--deg-engine",
        choices=["auto", "edger", "deseq2", "ttest"],
        default=None,
        help="组间 DEG 后端：auto 优先 rpy2 调用 edgeR，其次 DESeq2，否则 t-test+BH。任务描述里写 DESeq2/edgeR/t-test 也会被识别",
    )
    s.add_argument(
        "--marker-method",
        choices=["auto", "wilcoxon", "t-test", "mast"],
        default=None,
        help="探索性 cluster marker：Wilcoxon / t-test / MAST（MAST 需 R；不是组间结论）",
    )
    s.add_argument(
        "--deg-cross-validate",
        choices=["auto", "on", "off"],
        default=None,
        help="第二检验交叉验证基因列表。auto 默认开；任务里写「交叉验证」或「只用 Wilcoxon」可覆盖",
    )
    s.add_argument("--qc-method", choices=["mad", "percentile", "hybrid"], default=None)
    s.add_argument("--thread-id", default=None, help="LangGraph checkpoint thread_id；崩溃后用同一 id --resume")
    s.add_argument("--resume", action="store_true", help="从上次 checkpoint 续跑（默认同 --thread-id 或 .cache/last_thread_id）；会核对 run_manifest 的 scagent_version")
    s.add_argument("--from-checkpoint", action="store_true", dest="resume", help="同 --resume")
    s.add_argument("--force-resume", action="store_true", help="主版本不兼容时仍强制续跑")
    s.add_argument("--markers", default=None, help="自定义 marker CSV/JSON（可含 lineage 列）")
    s.add_argument("--report-lang", choices=["zh", "en", "both"], default="zh")
    s.add_argument("--selection", default=None, help="viewer 导出的 selection.json（框选细胞）")
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("view", help="Plotly 交互 UMAP：框选细胞")
    s.add_argument("--data", default=None, help="h5ad，默认 workspace/adata_processed.h5ad")
    s.add_argument("--serve", action="store_true", help="本地服务，框选后可直接提问")
    s.add_argument("--port", type=int, default=8765)
    s.add_argument("--open", action="store_true", help="用浏览器打开静态 viewer.html")
    s.set_defaults(func=cmd_view)

    s = sub.add_parser("ask", help="针对框选细胞提问（composition / marker in vs rest）")
    s.add_argument("query", nargs="?", default="", help="例如：分析我框选的这组细胞")
    s.add_argument("--selection", required=True, help="selection.json")
    s.add_argument("--data", default=None, help="h5ad")
    s.set_defaults(func=cmd_ask)

    s = sub.add_parser("confirm", help="确认 HITL 决策（线粒体阈值或 Leiden resolution）后继续")
    s.add_argument("kind", choices=["mt", "resolution"], help="决策点")
    s.add_argument("choice", help="lenient|recommended|strict 或 coarse|recommended|fine|0.4")
    s.add_argument("--execute", action="store_true")
    s.add_argument("--dry-run", action="store_true")
    s.add_argument("--thread-id", default=None)
    s.set_defaults(func=cmd_confirm)
    return p


def main(argv: list[str] | None = None) -> None:
    load_config()
    parser = build_parser()
    args = parser.parse_args(argv)
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
