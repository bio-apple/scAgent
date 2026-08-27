from __future__ import annotations

import json
import re

from agents.common import read_prompt, run_specialist
from agents.templates import (
    cluster_annotate_script,
    extract_locked_qc,
    qc_preprocess_script,
    splice_locked_qc,
)
from scagent.skills_loader import load_skill_text

PHASE_SKILLS = {
    "qc": ["anndata-data-structure", "scanpy-scrna-seq"],
    "downstream": [
        "scanpy-scrna-seq",
        "harmony-batch-correction",
        "scvi-tools-single-cell",
        "single-cell-annotation-guide",
        "celltypist-cell-annotation",
    ],
}


def _extract_code(text: str) -> str:
    m = re.search(r"```(?:python)?\n(.*?)```", text, re.DOTALL)
    if m:
        return m.group(1).strip()
    return text.strip()


def _skill_context(state: dict, phase: str) -> str:
    wanted = PHASE_SKILLS.get(phase) or []
    plan_skills = state.get("plan", {}).get("skills") or []
    names = [n for n in wanted if n in plan_skills] or wanted
    chunks = []
    for name in names:
        text = load_skill_text(name, include_references=False)
        chunks.append(text[:8000])
    return "\n\n".join(chunks)


def generate_code(state: dict, phase: str = "qc") -> str:
    meta = dict(state.get("metadata") or {})
    if state.get("markers_path"):
        meta["markers_path"] = state["markers_path"]
    qc = state.get("qc_strategy") or {}
    plan = state.get("plan") or {}
    if state.get("r_degraded") or plan.get("r_degraded"):
        return ""
    if phase == "qc":
        fallback = qc_preprocess_script(meta, qc)
        locked = extract_locked_qc(fallback) or ""
    else:
        fallback = cluster_annotate_script(meta, qc, plan)
        locked = ""
    review = state.get("review_qc" if phase == "qc" else "review_downstream") or state.get("review") or {}
    exe = state.get("execution_qc" if phase == "qc" else "execution_downstream") or state.get("execution") or {}
    arts = state.get("artifacts") or {}
    exec_feedback = {
        "ok": exe.get("ok"),
        "executed": exe.get("executed"),
        "returncode": exe.get("returncode"),
        "stderr_tail": (exe.get("stderr") or "")[-2000:],
        "stdout_tail": (exe.get("stdout") or "")[-1200:],
        "metrics": exe.get("metrics") or arts.get("metrics") or {},
        "issue_records": review.get("issue_records") or [],
        "issues": review.get("issues") or [],
    }
    llm = run_specialist(
        read_prompt("bio_coder"),
        (
            f"phase={phase}\n"
            f"metadata={json.dumps(meta, ensure_ascii=False)}\n"
            f"plan_route={plan.get('route')}\nintegrator={plan.get('integrator')}\n"
            f"qc={json.dumps({k: qc[k] for k in qc if k != 'rag_excerpt'}, ensure_ascii=False)}\n"
            f"reviewer_issues={json.dumps(review, ensure_ascii=False)}\n"
            f"execution_feedback={json.dumps(exec_feedback, ensure_ascii=False)}\n"
            f"skills:\n{_skill_context(state, phase)}\n"
            "输出完整可运行 Python。QC 阶段必须保留 LOCKED QC 块。"
            "注释阶段必须含 CellTypist + ≥2 阳性 + ≥1 阴性 marker；低置信度标 unknown/mixed。"
            "若 execution_feedback.ok 为 false：根据 stderr_tail 修复后再输出完整脚本。"
        ),
    )
    if not llm:
        return fallback
    code = _extract_code(llm)
    if phase == "qc":
        return splice_locked_qc(code, locked)
    return code
