from __future__ import annotations

import json
import re

from agents.common import read_prompt, run_specialist
from agents.templates import scanpy_script
from scagent.skills_loader import load_skill_text


def _extract_code(text: str) -> str:
    m = re.search(r"```(?:python)?\n(.*?)```", text, re.DOTALL)
    if m:
        return m.group(1).strip()
    return text.strip()


def generate_code(state: dict) -> str:
    meta = state.get("metadata") or {}
    qc = state.get("qc_strategy") or {}
    plan = state.get("plan") or {}
    review = state.get("review") or {}
    fallback = scanpy_script(meta, qc)
    skill_bits = []
    for name in (plan.get("skills") or [])[:3]:
        skill_bits.append(load_skill_text(name)[:4000])
    llm = run_specialist(
        read_prompt("bio_coder"),
        (
            f"metadata={json.dumps(meta, ensure_ascii=False)}\n"
            f"plan={json.dumps({k: plan[k] for k in plan if k != 'rag_excerpt'}, ensure_ascii=False)}\n"
            f"qc={json.dumps({k: qc[k] for k in qc if k != 'rag_excerpt'}, ensure_ascii=False)}\n"
            f"reviewer_issues={json.dumps(review, ensure_ascii=False)}\n"
            f"skills:\n{chr(10).join(skill_bits)}\n"
            "若有审查意见，请在代码中修复。"
        ),
    )
    if llm:
        return _extract_code(llm)
    return fallback
