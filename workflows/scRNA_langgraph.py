from __future__ import annotations

from langgraph.graph import END, START, StateGraph

from agents.annotation import build_annotation_plan
from agents.bio_coder import generate_code
from agents.planner import build_plan
from agents.qc_expert import build_qc_strategy
from agents.reviewer import review_state
from agents.writer import render_report
from sandbox.executor import write_and_maybe_run
from scagent.config import load_config, resolve_path
from scagent.inspect_data import inspect_data
from workflows.state import AgentState


def inspect_node(state: AgentState) -> dict:
    meta = inspect_data(
        state.get("data_path"),
        tissue=state.get("tissue"),
        extra={"task": state.get("user_query")},
    )
    if state.get("tissue"):
        meta["tissue"] = state["tissue"]
    return {"metadata": meta, "retry_count": state.get("retry_count") or 0}


def planner_node(state: AgentState) -> dict:
    plan = build_plan(state)
    return {"plan": plan, "skills_used": plan.get("skills") or []}


def qc_node(state: AgentState) -> dict:
    return {"qc_strategy": build_qc_strategy(state)}


def coder_node(state: AgentState) -> dict:
    return {"code": generate_code(state)}


def execute_node(state: AgentState) -> dict:
    cfg = load_config()
    workspace = resolve_path(cfg, "workspace")
    should = bool(state.get("execute_code"))
    result = write_and_maybe_run(
        state.get("code") or "",
        workspace=workspace,
        execute=should,
        timeout=int(cfg["analysis"]["timeout_seconds"]),
    )
    return {"execution": result}


def reviewer_node(state: AgentState) -> dict:
    return {"review": review_state(state)}


def annotation_node(state: AgentState) -> dict:
    return {"annotation_plan": build_annotation_plan(state)}


def writer_node(state: AgentState) -> dict:
    report = render_report(state)
    cfg = load_config()
    out_dir = resolve_path(cfg, "outputs")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "report.md"
    path.write_text(report, encoding="utf-8")
    return {"report": report}


def after_review(state: AgentState) -> str:
    cfg = load_config()
    review = state.get("review") or {}
    retries = int(state.get("retry_count") or 0)
    max_retries = int(cfg["analysis"]["max_review_retries"])
    if not review.get("passed") and retries < max_retries:
        return "retry"
    return "annotate"


def bump_retry(state: AgentState) -> dict:
    return {"retry_count": int(state.get("retry_count") or 0) + 1}


def build_graph():
    g = StateGraph(AgentState)
    g.add_node("inspect", inspect_node)
    g.add_node("planner", planner_node)
    g.add_node("qc_expert", qc_node)
    g.add_node("bio_coder", coder_node)
    g.add_node("execute", execute_node)
    g.add_node("reviewer", reviewer_node)
    g.add_node("retry_bump", bump_retry)
    g.add_node("annotation", annotation_node)
    g.add_node("writer", writer_node)

    g.add_edge(START, "inspect")
    g.add_edge("inspect", "planner")
    g.add_edge("planner", "qc_expert")
    g.add_edge("qc_expert", "bio_coder")
    g.add_edge("bio_coder", "execute")
    g.add_edge("execute", "reviewer")
    g.add_conditional_edges(
        "reviewer",
        after_review,
        {"retry": "retry_bump", "annotate": "annotation"},
    )
    g.add_edge("retry_bump", "bio_coder")
    g.add_edge("annotation", "writer")
    g.add_edge("writer", END)
    return g.compile()


def run_analysis(
    user_query: str,
    data_path: str | None = None,
    *,
    tissue: str = "default",
    language: str | None = None,
    execute_code: bool | None = None,
) -> AgentState:
    cfg = load_config()
    app = build_graph()
    state: AgentState = {
        "user_query": user_query,
        "data_path": data_path or "",
        "tissue": tissue,
        "language": language or cfg["analysis"]["language"],
        "execute_code": cfg["analysis"]["execute_code"] if execute_code is None else execute_code,
        "retry_count": 0,
    }
    return app.invoke(state)
