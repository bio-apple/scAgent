"""Plan-and-Solve multi-agent scRNA-seq graph: specialists + shared code audit/self-correct."""

from __future__ import annotations

from langgraph.graph import END, START, StateGraph

from agents.artifacts import skills_fingerprint
from scagent.compat import assert_resume_compatible, scagent_version
from agents.bio_interpret import build_interpretation_plan
from agents.cluster_deg import build_cluster_deg_plan
from agents.code_audit import generate_and_execute
from agents.memory import persist_memory
from agents.planner import build_plan
from agents.qc_expert import build_qc_strategy
from agents.reviewer import format_review_card, publication_review, review_state
from agents.writer import render_html, render_report, stage_report_figures, write_run_log
from scagent.config import load_config, resolve_path
from scagent.hitl import (
    build_mt_decision,
    build_resolution_decision,
    has_mt_choice,
    has_resolution_choice,
    need_hitl,
    pick_option,
    save_session,
    write_decision_card,
)
from scagent.inspect_data import inspect_data
from scagent.logutil import get_logger
from workflows.checkpointing import get_checkpointer, load_last_thread, remember_thread
from workflows.state import AgentState

log = get_logger("graph")


def _with_memory(state: AgentState, updates: dict) -> dict:
    merged = {**state, **updates}
    try:
        updates["analysis_memory"] = persist_memory(merged)
    except Exception as exc:
        log.warning("persist memory skipped: %s", exc)
    return updates


def _log(state: AgentState, msg: str) -> list[str]:
    logs = list(state.get("logs") or [])
    logs.append(msg)
    log.info("%s", msg)
    return logs


def inspect_node(state: AgentState) -> dict:
    logs = _log(state, "inspect: start")
    extra = {"task": state.get("user_query")}
    if state.get("batch_key"):
        extra["sample_key"] = state["batch_key"]
    meta = inspect_data(
        state.get("data_path"),
        tissue=state.get("tissue"),
        extra=extra,
    )
    if state.get("tissue"):
        meta["tissue"] = state["tissue"]
    if state.get("batch_key"):
        meta["sample_key"] = state["batch_key"]
        meta["need_batch_correction"] = True
    if state.get("markers_path"):
        meta["markers_path"] = state["markers_path"]
    logs = _log({**state, "logs": logs}, f"inspect: species={meta.get('species')} platform={meta.get('platform')} n_samples={meta.get('n_samples')}")
    return {
        "metadata": meta,
        "retry_count_qc": state.get("retry_count_qc") or 0,
        "retry_count_downstream": state.get("retry_count_downstream") or 0,
        "artifacts": {
            "skills_fingerprint": skills_fingerprint(),
            "scagent_version": scagent_version(),
            "phases": {},
        },
        "logs": logs,
        "status": "running",
    }


def planner_node(state: AgentState) -> dict:
    logs = _log(state, "planner: start")
    plan = build_plan(state)
    logs = _log({**state, "logs": logs}, f"planner: integrator={plan.get('integrator')} r_degraded={plan.get('r_degraded')}")
    return {
        "plan": plan,
        "skills_used": plan.get("skills") or [],
        "r_degraded": bool(plan.get("r_degraded")),
        "logs": logs,
    }


def hitl_mt_node(state: AgentState) -> dict:
    logs = _log(state, "hitl(mt): start")
    card = build_mt_decision(state)
    html = write_decision_card(card)
    save_session(state)
    wait = need_hitl(state) and not has_mt_choice(state)
    logs = _log(
        {**state, "logs": logs},
        f"hitl(mt): recommended={card.get('recommended')} wait={wait} card={html}",
    )
    return {"hitl_mt": card, "logs": logs}


def hitl_res_node(state: AgentState) -> dict:
    logs = _log(state, "hitl(resolution): start")
    card = build_resolution_decision(state)
    html = write_decision_card(card)
    updates: dict = {"hitl_resolution": card}
    if state.get("resolution") is None and has_resolution_choice(state):
        opt = pick_option(card, state.get("resolution_choice"))
        if opt.get("resolution") is not None:
            updates["resolution"] = float(opt["resolution"])
    save_session({**state, **updates})
    wait = need_hitl(state) and not has_resolution_choice({**state, **updates})
    logs = _log(
        {**state, "logs": logs},
        f"hitl(resolution): recommended={card.get('recommended')} wait={wait} card={html}",
    )
    updates["logs"] = logs
    return updates


def qc_node(state: AgentState) -> dict:
    logs = _log(state, "qc_preprocess: start")
    return {"qc_strategy": build_qc_strategy(state), "logs": logs, "phase": "qc"}


def code_audit_qc_node(state: AgentState) -> dict:
    logs = _log(state, "code_audit(qc): start")
    cfg = load_config()
    workspace = resolve_path(cfg, "workspace")
    updates = generate_and_execute(
        state,
        phase="qc",
        workspace=workspace,
        execute=bool(state.get("execute_code")),
        timeout=int(cfg["analysis"]["timeout_seconds"]),
    )
    logs = _log(
        {**state, "logs": logs},
        f"code_audit(qc): ok={updates.get('execution_qc', {}).get('ok')} executed={updates.get('execution_qc', {}).get('executed')} reused={updates.get('reused')}",
    )
    updates["logs"] = logs
    return _with_memory(state, updates)


def review_qc_node(state: AgentState) -> dict:
    logs = _log(state, "reviewer(qc): start")
    review = review_state({**state, "phase": "qc", "code": state.get("code_qc"), "execution": state.get("execution_qc")})
    logs = _log({**state, "logs": logs}, f"reviewer(qc): passed={review.get('passed')} issues={review.get('issues')}")
    return {"review_qc": review, "review": review, "logs": logs}


def bump_qc(state: AgentState) -> dict:
    n = int(state.get("retry_count_qc") or 0) + 1
    return {"retry_count_qc": n, "logs": _log(state, f"retry qc #{n}")}


def cluster_deg_node(state: AgentState) -> dict:
    logs = _log(state, "cluster_deg: start")
    plan = build_cluster_deg_plan(state)
    return {
        "cluster_deg_plan": plan,
        "annotation_plan": plan.get("annotation") or {},
        "logs": logs,
    }


def code_audit_down_node(state: AgentState) -> dict:
    logs = _log(state, "code_audit(cluster_deg): start")
    cfg = load_config()
    workspace = resolve_path(cfg, "workspace")
    updates = generate_and_execute(
        state,
        phase="downstream",
        workspace=workspace,
        execute=bool(state.get("execute_code")),
        timeout=int(cfg["analysis"]["timeout_seconds"]),
    )
    logs = _log(
        {**state, "logs": logs},
        f"code_audit(cluster_deg): ok={updates.get('execution_downstream', {}).get('ok')} reused={updates.get('reused')}",
    )
    updates["logs"] = logs
    return _with_memory(state, updates)


def review_down_node(state: AgentState) -> dict:
    logs = _log(state, "reviewer(cluster_deg): start")
    review = review_state(
        {**state, "phase": "downstream", "code": state.get("code_downstream"), "execution": state.get("execution_downstream")}
    )
    logs = _log({**state, "logs": logs}, f"reviewer(cluster_deg): passed={review.get('passed')}")
    return {"review_downstream": review, "review": review, "logs": logs}


def bump_down(state: AgentState) -> dict:
    n = int(state.get("retry_count_downstream") or 0) + 1
    return {"retry_count_downstream": n, "logs": _log(state, f"retry cluster_deg #{n}")}


def bio_interpret_node(state: AgentState) -> dict:
    logs = _log(state, "bio_interpret: start")
    plan = build_interpretation_plan(state)
    return {"interpretation_plan": plan, "phase": "interpret", "logs": logs}


def code_audit_interpret_node(state: AgentState) -> dict:
    logs = _log(state, "code_audit(interpret): start")
    cfg = load_config()
    workspace = resolve_path(cfg, "workspace")
    updates = generate_and_execute(
        state,
        phase="interpret",
        workspace=workspace,
        execute=bool(state.get("execute_code")),
        timeout=int(cfg["analysis"]["timeout_seconds"]),
    )
    logs = _log(
        {**state, "logs": logs},
        f"code_audit(interpret): ok={updates.get('execution_interpret', {}).get('ok')} reused={updates.get('reused')}",
    )
    updates["logs"] = logs
    return _with_memory(state, updates)


def review_pub_node(state: AgentState) -> dict:
    logs = _log(state, "reviewer(publication): start")
    card = publication_review(state)
    logs = _log(
        {**state, "logs": logs},
        f"reviewer(publication): score={card.get('score')}/{card.get('max_score')} verdict={card.get('verdict')}",
    )
    log.info("%s", format_review_card(card, state.get("report_lang") or "zh").strip())
    return {"review_publication": card, "logs": logs}


def writer_node(state: AgentState) -> dict:
    logs = _log(state, "writer: start")
    cfg = load_config()
    out_dir = resolve_path(cfg, "outputs")
    out_dir.mkdir(parents=True, exist_ok=True)
    mem = persist_memory(state, extra_dir=out_dir)
    arts = stage_report_figures(state.get("artifacts") or {}, out_dir)
    state_w = {**state, "analysis_memory": mem, "artifacts": arts}
    report = render_report(state_w)
    (out_dir / "report.md").write_text(report, encoding="utf-8")
    (out_dir / "report.html").write_text(render_html({**state_w, "report": report}), encoding="utf-8")
    write_run_log({**state_w, "report": report}, out_dir)
    from scagent.dual import export_dual
    from scagent.export_nb import export_analysis_notebook

    dual_path = export_dual(state_w, out_dir)
    nb_path = export_analysis_notebook(state_w, out_dir)
    logs = _log({**state, "logs": logs}, f"writer: dual={dual_path} notebook={nb_path}")
    from scagent.viewer import export_workspace_viewer

    viewer = None
    try:
        viewer = export_workspace_viewer(out_dir, state=state_w)
    except Exception as exc:
        log.warning("interactive viewer skipped: %s", exc)
    if viewer:
        logs = _log({**state, "logs": logs}, f"writer: viewer={viewer}")
    status = "completed"
    if state.get("r_degraded"):
        status = "r_degraded"
    elif need_hitl(state) and state.get("mode") != "annotate_only" and not has_mt_choice(state):
        status = "awaiting_mt_confirmation"
    elif need_hitl(state) and not has_resolution_choice(state):
        status = "awaiting_resolution_confirmation"
    elif state.get("interrupt_after_qc") and not state.get("auto_confirm") and not state.get("code_downstream"):
        status = "awaiting_qc_confirmation"
    elif state.get("mode") == "qc_only":
        status = "qc_only"
    elif state.get("execute_code") and (state.get("execution_qc") or {}).get("executed") and not (state.get("execution_qc") or {}).get("ok"):
        status = "qc_failed"
    logs = _log({**state, "logs": logs}, f"writer: wrote outputs/report.md status={status} notebook={nb_path.name}")
    out = {"report": report, "logs": logs, "status": status, "analysis_memory": mem, "notebook": str(nb_path), "artifacts": arts}
    if viewer:
        out["viewer"] = str(viewer)
    return out


def after_planner(state: AgentState) -> str:
    if state.get("r_degraded"):
        return "review_pub"
    if state.get("mode") == "annotate_only":
        return "hitl_res"
    return "hitl_mt"


def after_hitl_mt(state: AgentState) -> str:
    if need_hitl(state) and not has_mt_choice(state):
        return "review_pub"
    return "qc_expert"


def _run_failed(exe: dict, execute_code: bool) -> bool:
    if not execute_code:
        return False
    exe = exe or {}
    if exe.get("jail") in {"schema", "policy"}:
        return True
    return bool(exe.get("executed") and not exe.get("ok"))


def after_review_qc(state: AgentState) -> str:
    cfg = load_config()
    review = state.get("review_qc") or {}
    retries = int(state.get("retry_count_qc") or 0)
    max_retries = int(cfg["analysis"]["max_review_retries"])
    exe = state.get("execution_qc") or {}
    failed_run = _run_failed(exe, bool(state.get("execute_code")))
    not_passed = not review.get("passed")
    if (not_passed or failed_run) and retries < max_retries:
        return "retry"
    if state.get("mode") == "qc_only":
        return "review_pub"
    if failed_run or not_passed:
        return "review_pub"
    return "hitl_res"


def after_hitl_res(state: AgentState) -> str:
    if need_hitl(state) and not has_resolution_choice(state):
        return "review_pub"
    return "cluster_deg"


def after_review_down(state: AgentState) -> str:
    cfg = load_config()
    review = state.get("review_downstream") or {}
    retries = int(state.get("retry_count_downstream") or 0)
    max_retries = int(cfg["analysis"]["max_review_retries"])
    exe = state.get("execution_downstream") or {}
    failed_run = _run_failed(exe, bool(state.get("execute_code")))
    if (not review.get("passed") or failed_run) and retries < max_retries:
        return "retry"
    return "bio_interpret"


def build_graph(checkpointer=None):
    g = StateGraph(AgentState)
    g.add_node("inspect", inspect_node)
    g.add_node("planner", planner_node)
    g.add_node("hitl_mt", hitl_mt_node)
    g.add_node("qc_expert", qc_node)
    g.add_node("code_audit_qc", code_audit_qc_node)
    g.add_node("review_qc", review_qc_node)
    g.add_node("retry_qc", bump_qc)
    g.add_node("hitl_res", hitl_res_node)
    g.add_node("cluster_deg", cluster_deg_node)
    g.add_node("code_audit_down", code_audit_down_node)
    g.add_node("review_down", review_down_node)
    g.add_node("retry_down", bump_down)
    g.add_node("bio_interpret", bio_interpret_node)
    g.add_node("code_audit_interpret", code_audit_interpret_node)
    g.add_node("review_pub", review_pub_node)
    g.add_node("writer", writer_node)

    g.add_edge(START, "inspect")
    g.add_edge("inspect", "planner")
    g.add_conditional_edges(
        "planner",
        after_planner,
        {"hitl_mt": "hitl_mt", "hitl_res": "hitl_res", "review_pub": "review_pub"},
    )
    g.add_conditional_edges(
        "hitl_mt",
        after_hitl_mt,
        {"qc_expert": "qc_expert", "review_pub": "review_pub"},
    )
    g.add_edge("qc_expert", "code_audit_qc")
    g.add_edge("code_audit_qc", "review_qc")
    g.add_conditional_edges(
        "review_qc",
        after_review_qc,
        {"retry": "retry_qc", "hitl_res": "hitl_res", "review_pub": "review_pub"},
    )
    g.add_edge("retry_qc", "code_audit_qc")
    g.add_conditional_edges(
        "hitl_res",
        after_hitl_res,
        {"cluster_deg": "cluster_deg", "review_pub": "review_pub"},
    )
    g.add_edge("cluster_deg", "code_audit_down")
    g.add_edge("code_audit_down", "review_down")
    g.add_conditional_edges(
        "review_down",
        after_review_down,
        {"retry": "retry_down", "bio_interpret": "bio_interpret"},
    )
    g.add_edge("retry_down", "code_audit_down")
    g.add_edge("bio_interpret", "code_audit_interpret")
    g.add_edge("code_audit_interpret", "review_pub")
    g.add_edge("review_pub", "writer")
    g.add_edge("writer", END)
    cp = checkpointer if checkpointer is not None else get_checkpointer()
    return g.compile(checkpointer=cp)


def run_analysis(
    user_query: str,
    data_path: str | None = None,
    *,
    tissue: str = "default",
    language: str | None = None,
    execute_code: bool | None = None,
    mode: str = "full",
    interrupt_after_qc: bool = False,
    auto_confirm: bool = True,
    qc_choice: str | None = None,
    resolution_choice: str | None = None,
    resolution: float | None = None,
    batch_key: str | None = None,
    markers_path: str | None = None,
    report_lang: str = "zh",
    integrator: str | None = None,
    imputation: str | None = None,
    qc_method: str | None = None,
    remove_doublets: bool | None = None,
    doublet_methods: str | None = None,
    doublet_filter: str | None = None,
    ambient: str | None = None,
    condition_key: str | None = None,
    deg_engine: str | None = None,
    marker_method: str | None = None,
    deg_cross_validate: str | bool | None = None,
    thread_id: str | None = None,
    resume: bool = False,
    force_resume: bool = False,
    checkpointer=None,
    selection: dict | None = None,
) -> AgentState:
    from uuid import uuid4

    cfg = load_config()
    app = build_graph(checkpointer=checkpointer)
    if interrupt_after_qc and not auto_confirm:
        status_hint = "awaiting_mt_confirmation"
    else:
        status_hint = "running"
    tid = thread_id or (load_last_thread(cfg) if resume else None) or str(uuid4())
    lg_config = {"configurable": {"thread_id": tid}}
    remember_thread(tid, cfg)
    if resume:
        ws = resolve_path(cfg, "workspace")
        compat = assert_resume_compatible(ws / "run_manifest.json", force=force_resume)
        for w in compat.get("warnings") or []:
            log.warning("%s", w)
        updates: dict = {}
        if auto_confirm:
            updates["auto_confirm"] = True
        if qc_choice:
            updates["qc_choice"] = qc_choice
        if resolution_choice:
            updates["resolution_choice"] = resolution_choice
        if resolution is not None:
            updates["resolution"] = resolution
        if execute_code is not None:
            updates["execute_code"] = execute_code
        if mode:
            updates["mode"] = mode
        if updates:
            try:
                app.update_state(lg_config, updates)
            except Exception as exc:
                log.warning("checkpoint update_state skipped: %s", exc)
        out = app.invoke(None, lg_config)
        if out is not None:
            out["thread_id"] = tid
        return out
    state: AgentState = {
        "user_query": user_query,
        "data_path": data_path or "",
        "tissue": tissue,
        "language": language or cfg["analysis"]["language"],
        "execute_code": cfg["analysis"]["execute_code"] if execute_code is None else execute_code,
        "mode": mode,
        "interrupt_after_qc": interrupt_after_qc,
        "auto_confirm": auto_confirm,
        "qc_choice": qc_choice,
        "resolution_choice": resolution_choice,
        "resolution": resolution,
        "batch_key": batch_key,
        "markers_path": markers_path,
        "report_lang": report_lang,
        "integrator": integrator,
        "imputation": imputation,
        "qc_method": qc_method,
        "remove_doublets": remove_doublets,
        "doublet_methods": doublet_methods,
        "doublet_filter": doublet_filter,
        "ambient": ambient,
        "condition_key": condition_key,
        "deg_engine": deg_engine,
        "marker_method": marker_method,
        "deg_cross_validate": deg_cross_validate,
        "retry_count_qc": 0,
        "retry_count_downstream": 0,
        "logs": [],
        "status": status_hint,
        "thread_id": tid,
        "selection": selection or {},
    }
    out = app.invoke(state, lg_config)
    out["thread_id"] = tid
    return out
