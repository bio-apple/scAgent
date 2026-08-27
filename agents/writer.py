from __future__ import annotations

from datetime import date

from agents.common import read_prompt, run_specialist
from scagent.config import REPO_ROOT

TEMPLATE = REPO_ROOT / "report_templates" / "analysis_report.md"


def render_report(state: dict) -> str:
    tpl = TEMPLATE.read_text(encoding="utf-8") if TEMPLATE.exists() else "# 分析报告\n\n{body}\n"
    meta = state.get("metadata") or {}
    plan = state.get("plan") or {}
    qc = state.get("qc_strategy") or {}
    ann = state.get("annotation_plan") or {}
    review = state.get("review") or {}
    exe = state.get("execution") or {}
    figures = exe.get("figures") or []
    fig_md = "\n".join(f"- `{p}`" for p in figures) if figures else "- （本次未执行代码或未捕获图像）"

    body = tpl.format(
        date=date.today().isoformat(),
        query=state.get("user_query") or "",
        species=meta.get("species"),
        platform=meta.get("platform"),
        tissue=meta.get("tissue"),
        n_samples=meta.get("n_samples"),
        narrative=plan.get("narrative") or "",
        skills=", ".join(plan.get("skills") or []),
        route=" → ".join(plan.get("route") or []),
        qc_protocol=qc.get("protocol") or "",
        annotation=ann.get("instructions") or "",
        review_passed=review.get("passed"),
        review_issues="；".join(review.get("issues") or []) or "无",
        figures=fig_md,
        execution_ok=exe.get("ok"),
        stderr=(exe.get("stderr") or "")[:1500],
    )
    llm = run_specialist(
        read_prompt("writer"),
        f"请基于以下草稿写最终中文报告，不要编造未出现的结果。\n\n{body}",
    )
    return llm or body
