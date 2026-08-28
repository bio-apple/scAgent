from __future__ import annotations

import argparse
import json
import sys

from rag.ingest import ingest
from rag.retriever import clear_retrieve_cache, format_hits, retrieve, retrieve_fused
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
    cols = [c.strip() for c in args.collections.split(",") if c.strip()] if args.collections else None
    if cols:
        hits = retrieve(args.query, collections=cols, top_k=args.top_k)
    elif args.collection:
        hits = retrieve(args.query, collection=args.collection, top_k=args.top_k)
    else:
        hits = retrieve_fused(args.query, top_k=args.top_k)
    print(format_hits(hits))
    return 0


def cmd_skills(_args: argparse.Namespace) -> int:
    skills = list_skills()
    print(f"{len(skills)} bundled skills\n")
    for s in skills:
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
    from scagent.config import apply_performance_overrides
    from scagent.viewer import load_selection

    cfg = load_config()
    if getattr(args, "dask", False) or getattr(args, "gpu", False) or getattr(args, "rapids", False):
        apply_performance_overrides(
            cfg,
            dask=True if getattr(args, "dask", False) else None,
            gpu=True if getattr(args, "gpu", False) else None,
            rapids=True if getattr(args, "rapids", False) else None,
        )
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
            doublet_filter=args.doublet_filter,
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
    p = argparse.ArgumentParser(
        prog="scagent",
        description="LangGraph single-cell RNA-seq analysis agent. See README.en.md for full CLI reference.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("init", help="Interactive setup wizard (data path, tissue, task, resource limits)")
    s.add_argument(
        "--data",
        default=None,
        help="h5ad / loom / Cell Ranger outs / 10x mtx / Seurat .rds; comma-separated paths or sample parent dir",
    )
    s.add_argument("--tissue", default=None, help="Tissue type for QC profiles (e.g. pbmc, brain)")
    s.add_argument(
        "--task",
        choices=["annotate", "qc", "deg", "trajectory", "custom"],
        default=None,
        help="Preset task type",
    )
    s.add_argument("--query", default=None, help="Task description; default generated from --task")
    s.add_argument("--execute", action="store_true", help="Execute generated scripts after setup")
    s.add_argument("--dry-run", action="store_true", help="Generate scripts only (default)")
    s.add_argument("--yes", action="store_true", help="Accept defaults without prompts")
    s.add_argument("--run", action="store_true", help="Run scagent run immediately after writing config")
    s.add_argument(
        "--language",
        choices=["r_first", "python", "r"],
        default=None,
        help="r_first=R pipelines with Scanpy fallback; python=Scanpy only; r=export Rmd only",
    )
    s.add_argument("--report-lang", choices=["zh", "en", "both"], default=None, help="Report language")
    s.add_argument("--memory-mb", type=int, default=None, help="Sandbox memory limit (MB)")
    s.add_argument("--n-jobs", type=int, default=None, help="Scanpy/joblib parallelism; -1=all CPUs")
    s.add_argument("--timeout", type=int, default=None, help="Per-phase timeout (seconds)")
    s.add_argument("--batch-key", default=None, help="obs column name for batch/sample")
    s.add_argument("--condition-key", default=None, help="obs column for group comparison (pseudobulk DE)")
    s.add_argument("--interrupt", action="store_true", help="Enable HITL pauses (mito threshold, resolution)")
    s.add_argument("--no-write-config", action="store_true", help="Do not write config.local.yaml")
    s.set_defaults(func=cmd_init)

    s = sub.add_parser("ingest", help="Index knowledge/ (papers, methods, markers, best_practices, sops, upstream)")
    s.set_defaults(func=cmd_ingest)

    s = sub.add_parser(
        "update-kb",
        help="Fetch latest sc-best-practices (theislab) and rebuild RAG index",
    )
    s.add_argument(
        "--url",
        default=None,
        help="Git URL; default https://github.com/theislab/single-cell-best-practices.git",
    )
    s.add_argument(
        "--branch",
        default="main",
        help="Branch (default main; clone retries master on failure)",
    )
    s.set_defaults(func=cmd_update_kb)

    s = sub.add_parser("add-doc", help="Add lab SOP to local RAG (copy to knowledge/sops then ingest)")
    s.add_argument("path", help="File or directory (md / txt / pdf / ipynb)")
    s.add_argument("--name", default=None, help="Target filename or subdirectory under sops/")
    s.set_defaults(func=cmd_add_doc)

    s = sub.add_parser("retrieve", help="Query RAG collections")
    s.add_argument("query")
    s.add_argument("--collection", default=None, help="Single collection (default: fused across knowledge/)")
    s.add_argument("--collections", default=None, help="Comma-separated, e.g. papers,markers,best_practices")
    s.add_argument("--top-k", type=int, default=None, help="Number of chunks to return")
    s.set_defaults(func=cmd_retrieve)

    s = sub.add_parser("skills", help="List bundled SciAgent-style skills")
    s.set_defaults(func=cmd_skills)

    s = sub.add_parser("skill", help="Print one skill by name")
    s.add_argument("name")
    s.add_argument("--refs", action="store_true", help="Include reference links")
    s.set_defaults(func=cmd_show_skill)

    s = sub.add_parser("demo", help="Write 100-cell sparse demo h5ad")
    s.add_argument("--path", default=None, help="Output path (default tests/data/tiny_100cells.h5ad)")
    s.set_defaults(func=cmd_demo)

    s = sub.add_parser("memory", help="Print analysis provenance (steps and params, not chat history)")
    s.set_defaults(func=cmd_memory)

    s = sub.add_parser(
        "snapshots",
        help="List lazy AnnData snapshots for a thread (hardlinks / obs deltas, X not copied to RAM)",
    )
    s.add_argument("--thread-id", default=None, help="LangGraph thread_id (default last run)")
    s.set_defaults(func=cmd_snapshots)

    s = sub.add_parser(
        "branch",
        help="Fork a parameter experiment from a snapshot (shared h5ad inode, no matrix copy)",
    )
    s.add_argument("--from-thread", default=None, help="Source thread_id (default .cache/last_thread_id)")
    s.add_argument("--as", dest="as_name", required=True, help="New branch / thread name")
    s.add_argument("--step", default="qc", help="Fork point: qc or downstream")
    s.add_argument("--checkout", action="store_true", help="Hard-link snapshot into workspace/ for continued runs")
    s.set_defaults(func=cmd_branch)

    s = sub.add_parser("run", help="Planner → specialist agents → Code Audit → Reviewer → publication report")
    s.add_argument("query", help='Analysis task, e.g. "Standard PBMC QC and annotation"')
    s.add_argument(
        "--data",
        default="",
        help="h5ad / loom / Cell Ranger outs / 10x mtx / Seurat .rds; comma-separated = multi-sample concat",
    )
    s.add_argument("--tissue", default="default", help="Tissue profile for QC thresholds")
    s.add_argument(
        "--language",
        choices=["r_first", "python", "r"],
        default=None,
        help="r_first=R pipelines with Scanpy fallback; python=Scanpy only; r=export Rmd only",
    )
    s.add_argument("--execute", action="store_true", help="Execute generated scripts in workspace")
    s.add_argument("--dry-run", action="store_true", help="Generate scripts and report without execution")
    s.add_argument("--qc-only", action="store_true", help="Run QC phase only")
    s.add_argument(
        "--annotate-only",
        action="store_true",
        help="Skip QC; requires existing workspace/adata_qc.h5ad",
    )
    s.add_argument(
        "--interrupt",
        action="store_true",
        help="Pause at mitochondrial threshold and Leiden resolution for human confirmation",
    )
    s.add_argument("--yes", action="store_true", help="Skip human confirmation (use recommended presets)")
    s.add_argument("--resolution", type=float, default=None, help="Fixed Leiden resolution (skip HITL resolution step)")
    s.add_argument("--batch-key", default=None, help="obs column for batch integration")
    s.add_argument(
        "--integrator",
        choices=["auto", "none", "harmony", "scvi", "cca", "scanorama", "bbknn"],
        default=None,
        help="Batch correction: auto=Harmony or scVI after inspect; cca/scanorama=Scanorama; bbknn=neighbor graph",
    )
    s.add_argument(
        "--impute",
        choices=["none", "magic", "alra"],
        default=None,
        help="Dropout imputation (default from config.modules.imputation)",
    )
    s.add_argument(
        "--ambient",
        choices=["auto", "none", "soupx", "decontx"],
        default=None,
        help="Ambient RNA removal (default from config.modules.ambient)",
    )
    s.add_argument(
        "--remove-doublets",
        action="store_true",
        help="Filter doublets per --doublet-filter (default high_conf = high-confidence only)",
    )
    s.add_argument(
        "--doublet-filter",
        choices=["high_conf", "all"],
        default=None,
        help="high_conf=conservative (high-confidence only); all=strict (high + low confidence)",
    )
    s.add_argument(
        "--doublet-methods",
        choices=["auto", "scrublet", "both"],
        default=None,
        help="auto=Scrublet+scDblFinder on multi-sample/complex tissue; both=force cross-validation",
    )
    s.add_argument(
        "--condition-key",
        default=None,
        help="obs column for group comparison; enables sample-level pseudobulk DE",
    )
    s.add_argument(
        "--deg-engine",
        choices=["auto", "edger", "deseq2", "ttest"],
        default=None,
        help="Group DEG backend: auto prefers rpy2 edgeR, then DESeq2, then t-test+BH",
    )
    s.add_argument(
        "--marker-method",
        choices=["auto", "wilcoxon", "t-test", "mast"],
        default=None,
        help="Exploratory cluster markers (Wilcoxon / t-test / MAST); not for group-level conclusions",
    )
    s.add_argument(
        "--deg-cross-validate",
        choices=["auto", "on", "off"],
        default=None,
        help="Second statistical test to cross-validate DEG gene lists",
    )
    s.add_argument(
        "--qc-method",
        choices=["mad", "percentile", "hybrid"],
        default=None,
        help="Dynamic QC threshold method",
    )
    s.add_argument(
        "--thread-id",
        default=None,
        help="LangGraph checkpoint thread_id; reuse with --resume after crash",
    )
    s.add_argument(
        "--resume",
        action="store_true",
        help="Resume from last checkpoint (same --thread-id or .cache/last_thread_id); checks scagent_version",
    )
    s.add_argument("--from-checkpoint", action="store_true", dest="resume", help="Alias for --resume")
    s.add_argument(
        "--force-resume",
        action="store_true",
        help="Force resume even when major scagent_version mismatch",
    )
    s.add_argument("--markers", default=None, help="Custom marker CSV/JSON (optional lineage column)")
    s.add_argument("--report-lang", choices=["zh", "en", "both"], default="zh", help="Publication report language")
    s.add_argument("--selection", default=None, help="selection.json exported from viewer (selected cells)")
    s.add_argument("--dask", action="store_true", help="Enable experimental Dask/out-of-core path (large h5ad)")
    s.add_argument("--gpu", action="store_true", help="Enable GPU for scVI training when CUDA is available")
    s.add_argument(
        "--rapids",
        action="store_true",
        help="Use RAPIDS for neighbors/UMAP (requires rapids-singlecell; implies --gpu)",
    )
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("view", help="Interactive Plotly UMAP viewer with cell selection")
    s.add_argument("--data", default=None, help="h5ad path (default workspace/adata_processed.h5ad)")
    s.add_argument("--serve", action="store_true", help="Start local server for live Q&A after selection")
    s.add_argument("--port", type=int, default=8765, help="Port for --serve")
    s.add_argument("--open", action="store_true", help="Open static outputs/viewer.html in browser")
    s.set_defaults(func=cmd_view)

    s = sub.add_parser("ask", help="Ask about lasso/box-selected cells (composition / markers vs rest)")
    s.add_argument("query", nargs="?", default="", help='e.g. "Analyze my selected cells"')
    s.add_argument("--selection", required=True, help="Path to selection.json from viewer")
    s.add_argument("--data", default=None, help="h5ad path")
    s.set_defaults(func=cmd_ask)

    s = sub.add_parser("confirm", help="Confirm HITL decision (mito threshold or Leiden resolution) and continue")
    s.add_argument("kind", choices=["mt", "resolution"], help="Decision point: mt or resolution")
    s.add_argument(
        "choice",
        help="mt: lenient|recommended|strict; resolution: coarse|recommended|fine|<float>",
    )
    s.add_argument("--execute", action="store_true", help="Execute scripts after confirm")
    s.add_argument("--dry-run", action="store_true", help="Generate only after confirm")
    s.add_argument("--thread-id", default=None, help="LangGraph thread_id")
    s.set_defaults(func=cmd_confirm)
    return p


def main(argv: list[str] | None = None) -> None:
    load_config()
    parser = build_parser()
    args = parser.parse_args(argv)
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
