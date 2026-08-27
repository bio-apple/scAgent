from __future__ import annotations

import hashlib

from langgraph.graph import END, START, StateGraph

from agents.annotation import build_annotation_plan
from agents.artifacts import collect_workspace, merge_artifacts, skills_fingerprint
from agents.bio_coder import generate_code
from agents.memory import persist_memory
from agents.planner import build_plan
from agents.qc_expert import build_qc_strategy
from agents.reviewer import format_review_card, publication_review, review_state
from agents.writer import render_html, render_report, write_run_log
from sandbox.executor import write_and_maybe_run
from scagent.config import load_config, resolve_path
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
        "artifacts": {"skills_fingerprint": skills_fingerprint(), "phases": {}},
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


def qc_node(state: AgentState) -> dict:
    logs = _log(state, "qc_expert: start")
    return {"qc_strategy": build_qc_strategy(state), "logs": logs, "phase": "qc"}


def coder_qc_node(state: AgentState) -> dict:
    logs = _log(state, "bio_coder(qc): start")
    code = generate_code(state, phase="qc")
    return {"code_qc": code, "code": code, "phase": "qc", "logs": logs}


def _code_fp(code: str | None) -> str:
    return hashlib.sha256((code or "").encode("utf-8")).hexdigest()[:16]


def _reuse_execution(prev: dict | None, code: str | None, want_exec: bool) -> dict | None:
    """Reuse a successful checkpointed run unless code changed or execute was skipped before."""
    prev = prev or {}
    if not prev.get("ok"):
        return None
    if prev.get("code_fp") != _code_fp(code):
        return None
    if want_exec and not prev.get("executed"):
        return None
    return prev


def execute_qc_node(state: AgentState) -> dict:
    cfg = load_config()
    workspace = resolve_path(cfg, "workspace")
    want_exec = bool(state.get("execute_code"))
    reused = _reuse_execution(state.get("execution_qc"), state.get("code_qc"), want_exec)
    if reused is not None:
        logs = _log(state, "execute(qc): skip (checkpoint ok)")
        return _with_memory(state, {"execution_qc": reused, "execution": reused, "logs": logs})
    logs = _log(state, "execute(qc): start")
    result = write_and_maybe_run(
        state.get("code_qc") or "",
        workspace=workspace,
        execute=want_exec,
        timeout=int(cfg["analysis"]["timeout_seconds"]),
        filename="qc_preprocess.py",
        extra_manifest={"phase": "qc", "data_path": state.get("data_path")},
    )
    result = dict(result)
    result["code_fp"] = _code_fp(state.get("code_qc"))
    art = merge_artifacts(state.get("artifacts"), collect_workspace(workspace, "qc", result))
    logs = _log({**state, "logs": logs}, f"execute(qc): ok={result.get('ok')} executed={result.get('executed')}")
    return _with_memory(
        state,
        {"execution_qc": result, "execution": result, "artifacts": art, "logs": logs},
    )


def review_qc_node(state: AgentState) -> dict:
    logs = _log(state, "reviewer(qc): start")
    review = review_state({**state, "phase": "qc", "code": state.get("code_qc"), "execution": state.get("execution_qc")})
    logs = _log({**state, "logs": logs}, f"reviewer(qc): passed={review.get('passed')} issues={review.get('issues')}")
    return {"review_qc": review, "review": review, "logs": logs}


def bump_qc(state: AgentState) -> dict:
    n = int(state.get("retry_count_qc") or 0) + 1
    return {"retry_count_qc": n, "logs": _log(state, f"retry qc #{n}")}


def annotation_node(state: AgentState) -> dict:
    logs = _log(state, "annotation: start")
    plan = build_annotation_plan(state)
    return {"annotation_plan": plan, "logs": logs}


def coder_down_node(state: AgentState) -> dict:
    logs = _log(state, "bio_coder(downstream): start")
    ann = state.get("annotation_plan") or build_annotation_plan(state)
    code = generate_code(state, phase="downstream") or ann.get("code") or ""
    return {
        "annotation_plan": ann,
        "code_downstream": code,
        "code": code,
        "phase": "downstream",
        "logs": logs,
    }


def execute_down_node(state: AgentState) -> dict:
    cfg = load_config()
    workspace = resolve_path(cfg, "workspace")
    want_exec = bool(state.get("execute_code"))
    reused = _reuse_execution(state.get("execution_downstream"), state.get("code_downstream"), want_exec)
    if reused is not None:
        logs = _log(state, "execute(downstream): skip (checkpoint ok)")
        return _with_memory(
            state, {"execution_downstream": reused, "execution": reused, "logs": logs}
        )
    logs = _log(state, "execute(downstream): start")
    result = write_and_maybe_run(
        state.get("code_downstream") or "",
        workspace=workspace,
        execute=want_exec,
        timeout=int(cfg["analysis"]["timeout_seconds"]),
        filename="cluster_annotate.py",
        extra_manifest={"phase": "downstream", "data_path": state.get("data_path")},
    )
    result = dict(result)
    result["code_fp"] = _code_fp(state.get("code_downstream"))
    art = merge_artifacts(state.get("artifacts"), collect_workspace(workspace, "downstream", result))
    logs = _log({**state, "logs": logs}, f"execute(downstream): ok={result.get('ok')}")
    return _with_memory(
        state,
        {"execution_downstream": result, "execution": result, "artifacts": art, "logs": logs},
    )


def review_down_node(state: AgentState) -> dict:
    logs = _log(state, "reviewer(downstream): start")
    review = review_state(
        {**state, "phase": "downstream", "code": state.get("code_downstream"), "execution": state.get("execution_downstream")}
    )
    logs = _log({**state, "logs": logs}, f"reviewer(downstream): passed={review.get('passed')}")
    return {"review_downstream": review, "review": review, "logs": logs}


def bump_down(state: AgentState) -> dict:
    n = int(state.get("retry_count_downstream") or 0) + 1
    return {"retry_count_downstream": n, "logs": _log(state, f"retry downstream #{n}")}


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
    state_w = {**state, "analysis_memory": mem}
    report = render_report(state_w)
    (out_dir / "report.md").write_text(report, encoding="utf-8")
    (out_dir / "report.html").write_text(render_html({**state_w, "report": report}), encoding="utf-8")
    write_run_log({**state_w, "report": report}, out_dir)
    status = "completed"
    if state.get("r_degraded"):
        status = "r_degraded"
    elif state.get("interrupt_after_qc") and not state.get("auto_confirm") and not state.get("code_downstream"):
        status = "awaiting_qc_confirmation"
    elif state.get("mode") == "qc_only":
        status = "qc_only"
    elif state.get("execute_code") and (state.get("execution_qc") or {}).get("executed") and not (state.get("execution_qc") or {}).get("ok"):
        status = "qc_failed"
    logs = _log({**state, "logs": logs}, f"writer: wrote outputs/report.md status={status}")
    return {"report": report, "logs": logs, "status": status, "analysis_memory": mem}


def after_planner(state: AgentState) -> str:
    if state.get("r_degraded"):
        return "review_pub"
    if state.get("mode") == "annotate_only":
        return "annotation"
    return "qc_expert"


def after_review_qc(state: AgentState) -> str:
    cfg = load_config()
    review = state.get("review_qc") or {}
    retries = int(state.get("retry_count_qc") or 0)
    max_retries = int(cfg["analysis"]["max_review_retries"])
    exe = state.get("execution_qc") or {}
    failed_run = bool(state.get("execute_code") and exe.get("executed") and not exe.get("ok"))
    not_passed = not review.get("passed")
    if (not_passed or failed_run) and retries < max_retries:
        return "retry"
    if state.get("mode") == "qc_only":
        return "review_pub"
    if state.get("interrupt_after_qc") and not state.get("auto_confirm"):
        return "review_pub"
    if failed_run or not_passed:
        return "review_pub"
    return "annotation"


def after_review_down(state: AgentState) -> str:
    cfg = load_config()
    review = state.get("review_downstream") or {}
    retries = int(state.get("retry_count_downstream") or 0)
    max_retries = int(cfg["analysis"]["max_review_retries"])
    exe = state.get("execution_downstream") or {}
    failed_run = bool(state.get("execute_code") and exe.get("executed") and not exe.get("ok"))
    if (not review.get("passed") or failed_run) and retries < max_retries:
        return "retry"
    return "review_pub"


def build_graph(checkpointer=None):
    g = StateGraph(AgentState)
    g.add_node("inspect", inspect_node)
    g.add_node("planner", planner_node)
    g.add_node("qc_expert", qc_node)
    g.add_node("coder_qc", coder_qc_node)
    g.add_node("execute_qc", execute_qc_node)
    g.add_node("review_qc", review_qc_node)
    g.add_node("retry_qc", bump_qc)
    g.add_node("annotation", annotation_node)
    g.add_node("coder_down", coder_down_node)
    g.add_node("execute_down", execute_down_node)
    g.add_node("review_down", review_down_node)
    g.add_node("retry_down", bump_down)
    g.add_node("review_pub", review_pub_node)
    g.add_node("writer", writer_node)

    g.add_edge(START, "inspect")
    g.add_edge("inspect", "planner")
    g.add_conditional_edges(
        "planner",
        after_planner,
        {"qc_expert": "qc_expert", "annotation": "annotation", "review_pub": "review_pub"},
    )
    g.add_edge("qc_expert", "coder_qc")
    g.add_edge("coder_qc", "execute_qc")
    g.add_edge("execute_qc", "review_qc")
    g.add_conditional_edges(
        "review_qc",
        after_review_qc,
        {"retry": "retry_qc", "annotation": "annotation", "review_pub": "review_pub"},
    )
    g.add_edge("retry_qc", "coder_qc")
    g.add_edge("annotation", "coder_down")
    g.add_edge("coder_down", "execute_down")
    g.add_edge("execute_down", "review_down")
    g.add_conditional_edges(
        "review_down",
        after_review_down,
        {"retry": "retry_down", "review_pub": "review_pub"},
    )
    g.add_edge("retry_down", "coder_down")
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
    resolution: float | None = None,
    batch_key: str | None = None,
    markers_path: str | None = None,
    report_lang: str = "zh",
    integrator: str | None = None,
    imputation: str | None = None,
    qc_method: str | None = None,
    remove_doublets: bool | None = None,
    ambient: str | None = None,
    condition_key: str | None = None,
    thread_id: str | None = None,
    resume: bool = False,
    checkpointer=None,
) -> AgentState:
    from uuid import uuid4

    cfg = load_config()
    app = build_graph(checkpointer=checkpointer)
    if interrupt_after_qc and not auto_confirm:
        status_hint = "awaiting_qc_confirmation"
    else:
        status_hint = "running"
    tid = thread_id or (load_last_thread(cfg) if resume else None) or str(uuid4())
    lg_config = {"configurable": {"thread_id": tid}}
    remember_thread(tid, cfg)
    if resume:
        updates: dict = {}
        if auto_confirm:
            updates["auto_confirm"] = True
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
        if interrupt_after_qc and not auto_confirm and (out or {}).get("review_qc", {}).get("passed"):
            out["status"] = "awaiting_qc_confirmation"
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
        "resolution": resolution,
        "batch_key": batch_key,
        "markers_path": markers_path,
        "report_lang": report_lang,
        "integrator": integrator,
        "imputation": imputation,
        "qc_method": qc_method,
        "remove_doublets": remove_doublets,
        "ambient": ambient,
        "condition_key": condition_key,
        "retry_count_qc": 0,
        "retry_count_downstream": 0,
        "logs": [],
        "status": status_hint,
        "thread_id": tid,
    }
    out = app.invoke(state, lg_config)
    if interrupt_after_qc and not auto_confirm and out.get("review_qc", {}).get("passed"):
        out["status"] = "awaiting_qc_confirmation"
    out["thread_id"] = tid
    return out
