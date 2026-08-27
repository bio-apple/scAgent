from __future__ import annotations

from typing import Any, TypedDict


class AgentState(TypedDict, total=False):
    user_query: str
    data_path: str
    tissue: str
    language: str
    execute_code: bool
    metadata: dict[str, Any]
    plan: dict[str, Any]
    qc_strategy: dict[str, Any]
    code: str
    execution: dict[str, Any]
    annotation_plan: dict[str, Any]
    review: dict[str, Any]
    retry_count: int
    report: str
    skills_used: list[str]
