from __future__ import annotations

import json
import re

from agents.code_schema import validate_script
from agents.common import read_prompt, run_specialist
from agents.templates import (
    cluster_annotate_script,
    extract_locked_qc,
    qc_preprocess_script,
    splice_locked_qc,
)
from scagent.logutil import get_logger
from scagent.deg_methods import force_pseudobulk_de
from scagent.skills_loader import load_skill_text

log = get_logger("bio_coder")

PHASE_SKILLS = {
    "qc": ["anndata-data-structure", "scanpy-scrna-seq"],
    "downstream": [
        "scanpy-scrna-seq",
        "harmony-batch-correction",
        "scvi-tools-single-cell",
        "single-cell-annotation-guide",
        "single-cell-annotation",
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


def _exec_feedback(state: dict, phase: str) -> dict:
    review = state.get("review_qc" if phase == "qc" else "review_downstream") or state.get("review") or {}
    exe = state.get("execution_qc" if phase == "qc" else "execution_downstream") or state.get("execution") or {}
    arts = state.get("artifacts") or {}
    return {
        "ok": exe.get("ok"),
        "executed": exe.get("executed"),
        "returncode": exe.get("returncode"),
        "stderr_tail": (exe.get("stderr") or "")[-2000:],
        "stdout_tail": (exe.get("stdout") or "")[-1200:],
        "metrics": exe.get("metrics") or arts.get("metrics") or {},
        "issue_records": review.get("issue_records") or [],
        "issues": review.get("issues") or [],
        "schema": exe.get("schema") or {},
        "jail": exe.get("jail"),
    }


def _deg_hard_rules(state: dict, phase: str) -> str:
    if phase != "downstream":
        return ""
    meta = state.get("metadata") or {}
    plan = state.get("plan") or {}
    if not force_pseudobulk_de(meta, plan):
        return ""
    ck = meta.get("condition_key") or plan.get("condition_key") or "condition"
    return (
        f"\n【硬约束】condition_key={ck!r} 且 n_replicates≥2："
        "组间差异表达 MUST 调用 pseudobulk_de（sample×cell_type + DESeq2/edgeR + FDR）。"
        f"禁止 sc.tl.rank_genes_groups(groupby={ck!r}) 或任何 cell-level Wilcoxon/MAST 作为组间结论。"
        "rank_genes/leiden 探索性 cluster marker 仍可用 Wilcoxon。"
    )


def _coder_task(state: dict, phase: str, *, extra: str = "") -> str:
    meta = dict(state.get("metadata") or {})
    if state.get("markers_path"):
        meta["markers_path"] = state["markers_path"]
    qc = state.get("qc_strategy") or {}
    plan = dict(state.get("plan") or {})
    if state.get("resolution") is not None:
        plan["resolution"] = state["resolution"]
    dag = plan.get("dag") or {}
    return (
        f"phase={phase}\n"
        f"metadata={json.dumps(meta, ensure_ascii=False)}\n"
        f"plan_route={plan.get('route')}\n"
        f"dag={json.dumps(dag, ensure_ascii=False)}\n"
        f"integrator={plan.get('integrator')}\n"
        f"qc={json.dumps({k: qc[k] for k in qc if k != 'rag_excerpt'}, ensure_ascii=False)}\n"
        f"reviewer_issues={json.dumps(state.get('review_qc' if phase == 'qc' else 'review_downstream') or state.get('review') or {}, ensure_ascii=False)}\n"
        f"execution_feedback={json.dumps(_exec_feedback(state, phase), ensure_ascii=False)}\n"
        f"skills:\n{_skill_context(state, phase)}\n"
        + (f"tool_route={json.dumps((state.get('plan') or {}).get('tool_route') or {}, ensure_ascii=False)}\n")
        + "Always use R ecosystem first. Only invoke Python when R lacks the required functionality.\n"
        "输出完整可运行 Python。QC 阶段必须保留 LOCKED QC 块。"
        "注释阶段必须含 CellTypist + scANVI 集成（ensemble_cell_annotation → obs['scagent_annotation']）+ ≥2 阳性 + ≥1 阴性 marker，并用 fuse_annotation 融合；禁止只调用 Azimuth。"
        "必须遵守 DAG：PCA/neighbors/UMAP/Leiden 之后才能 rank_genes_groups、pseudobulk_de 或 DPT/PAGA。"
        + _deg_hard_rules(state, phase)
        + "若 execution_feedback.ok 为 false：根据 stderr_tail 修复语法或参数后再输出完整脚本。"
        + extra
    )


def _finalize(code: str, phase: str, locked: str) -> str:
    if phase == "qc":
        return splice_locked_qc(code, locked)
    return code


def generate_code(state: dict, phase: str = "qc") -> str:
    meta = dict(state.get("metadata") or {})
    if state.get("markers_path"):
        meta["markers_path"] = state["markers_path"]
    qc = state.get("qc_strategy") or {}
    plan = dict(state.get("plan") or {})
    if state.get("resolution") is not None:
        plan["resolution"] = state["resolution"]
    if state.get("r_degraded") or plan.get("r_degraded"):
        return ""
    if phase == "qc":
        fallback = qc_preprocess_script(meta, qc)
        locked = extract_locked_qc(fallback) or ""
    else:
        fallback = cluster_annotate_script(meta, qc, plan)
        locked = ""
    lang = str(state.get("language") or plan.get("language") or "python")
    llm = run_specialist(read_prompt("bio_coder"), _coder_task(state, phase))
    if not llm:
        return fallback
    code = _finalize(_extract_code(llm), phase, locked)
    report = validate_script(code, phase=phase, language=lang, metadata=meta, plan=plan)
    if report.get("ok"):
        return code
    log.warning("schema rejected generated %s code: %s", phase, report.get("issues"))
    extra = (
        "\n上一版未通过 schema/DAG 校验: "
        + json.dumps(report.get("issues"), ensure_ascii=False)
        + "\n请修复调用顺序或语法后重新输出完整脚本，不要只输出 diff。"
    )
    llm2 = run_specialist(read_prompt("bio_coder"), _coder_task(state, phase, extra=extra))
    if llm2:
        code2 = _finalize(_extract_code(llm2), phase, locked)
        report2 = validate_script(code2, phase=phase, language=lang, metadata=meta, plan=plan)
        if report2.get("ok"):
            return code2
        log.warning("self-correct still failed schema: %s", report2.get("issues"))
    return fallback
