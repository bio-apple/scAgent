"""Code Audit & Execution Agent: specialist instructions → sandbox script → schema gate → self-correct."""

from __future__ import annotations

import hashlib

from agents.bio_coder import generate_code
from agents.templates import annotate_deg_script, cluster_only_script, interpret_pathways_script
from sandbox.executor import write_and_maybe_run
from scagent.logutil import get_logger

log = get_logger("code_audit")

PHASE_IO = {
    "qc": ("qc_preprocess.py", "code_qc", "execution_qc"),
    "downstream": ("cluster_annotate.py", "code_downstream", "execution_downstream"),
    "interpret": ("interpret_pathways.py", "code_interpret", "execution_interpret"),
}

DOWNSTREAM_SPLIT = (
    ("cluster_only.py", "code_cluster", "execution_cluster"),
    ("annotate_deg.py", "code_annotate", "execution_annotate"),
)


def _downstream_plan(state: dict) -> dict:
    plan = dict(state.get("plan") or state.get("cluster_deg_plan") or {})
    if state.get("resolution") is not None:
        plan["resolution"] = state["resolution"]
    return plan


def _code_fp(code: str | None) -> str:
    return hashlib.sha256((code or "").encode("utf-8")).hexdigest()[:16]


def _reuse(prev: dict | None, code: str | None, want_exec: bool) -> dict | None:
    prev = prev or {}
    if not prev.get("ok"):
        return None
    if prev.get("code_fp") != _code_fp(code):
        return None
    if want_exec and not prev.get("executed"):
        return None
    return prev


def _script_for(state: dict, phase: str) -> str:
    if phase == "interpret":
        return interpret_pathways_script(
            state.get("metadata") or {},
            state.get("plan") or {},
            state.get("interpretation_plan") or {},
        )
    if phase == "qc":
        return generate_code(state, phase="qc")
    if phase == "cluster":
        meta = state.get("metadata") or {}
        qc = state.get("qc_strategy") or {}
        plan = _downstream_plan(state)
        return cluster_only_script(meta, qc, plan)
    if phase == "annotate":
        meta = state.get("metadata") or {}
        qc = state.get("qc_strategy") or {}
        plan = _downstream_plan(state)
        return annotate_deg_script(meta, qc, plan)
    ann = state.get("annotation_plan") or {}
    return generate_code(state, phase="downstream") or ann.get("code") or ""


def generate_and_execute(
    state: dict,
    *,
    phase: str,
    workspace,
    execute: bool,
    timeout: int,
) -> dict:
    """One agent turn: emit code, jail/schema, run, attach artifacts. Self-correct lives in generate_code + graph retry."""
    from agents.artifacts import collect_workspace, merge_artifacts

    if phase == "downstream":
        return _generate_downstream_split(state, workspace=workspace, execute=execute, timeout=timeout)

    filename, code_key, exec_key = PHASE_IO[phase]
    code = _script_for(state, phase)
    reused = _reuse(state.get(exec_key), code, execute)
    if reused is not None:
        log.info("code_audit(%s): reuse checkpoint", phase)
        return {
            code_key: code,
            "code": code,
            exec_key: reused,
            "execution": reused,
            "phase": phase,
            "reused": True,
        }
    result = write_and_maybe_run(
        code or "",
        workspace=workspace,
        execute=execute,
        timeout=timeout,
        filename=filename,
        extra_manifest={"phase": phase, "data_path": state.get("data_path"), "thread_id": state.get("thread_id")},
    )
    result = dict(result)
    result["code_fp"] = _code_fp(code)
    art = merge_artifacts(state.get("artifacts"), collect_workspace(workspace, phase, result))
    log.info("code_audit(%s): ok=%s executed=%s jail=%s", phase, result.get("ok"), result.get("executed"), result.get("jail"))
    return {
        code_key: code,
        "code": code,
        exec_key: result,
        "execution": result,
        "artifacts": art,
        "phase": phase,
        "reused": False,
    }


def _generate_downstream_split(state: dict, *, workspace, execute: bool, timeout: int) -> dict:
    from agents.artifacts import collect_workspace, merge_artifacts
    from agents.templates import cluster_annotate_script

    meta = state.get("metadata") or {}
    qc = state.get("qc_strategy") or {}
    plan = _downstream_plan(state)
    combined = cluster_annotate_script(meta, qc, plan, phase="full")
    cluster_code = cluster_only_script(meta, qc, plan)
    annotate_code = annotate_deg_script(meta, qc, plan)
    updates: dict = {"phase": "downstream", "reused": False}
    last_result = None
    for filename, code_key, exec_key in DOWNSTREAM_SPLIT:
        code = cluster_code if filename == "cluster_only.py" else annotate_code
        updates[code_key] = code
        reused = _reuse(state.get(exec_key), code, execute)
        if reused is not None:
            updates[exec_key] = reused
            last_result = reused
            continue
        result = write_and_maybe_run(
            code,
            workspace=workspace,
            execute=execute,
            timeout=timeout,
            filename=filename,
            extra_manifest={"phase": filename.replace(".py", ""), "data_path": state.get("data_path"), "thread_id": state.get("thread_id")},
        )
        result = dict(result)
        result["code_fp"] = _code_fp(code)
        updates[exec_key] = result
        last_result = result
    updates["code_downstream"] = combined
    updates["code"] = combined
    (workspace / "cluster_annotate.py").write_text(combined, encoding="utf-8")
    if last_result is not None:
        updates["execution_downstream"] = last_result
        updates["execution"] = last_result
        updates["artifacts"] = merge_artifacts(state.get("artifacts"), collect_workspace(workspace, "downstream", last_result))
    return updates
