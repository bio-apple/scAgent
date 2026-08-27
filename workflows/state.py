from __future__ import annotations

from typing import Any, TypedDict


class AgentState(TypedDict, total=False):
    user_query: str
    data_path: str
    tissue: str
    language: str
    execute_code: bool
    mode: str  # full | qc_only | annotate_only
    interrupt_after_qc: bool
    auto_confirm: bool
    qc_choice: str | None
    resolution_choice: str | None
    hitl_mt: dict[str, Any]
    hitl_resolution: dict[str, Any]
    resolution: float | None
    batch_key: str | None
    markers_path: str | None
    integrator: str | None
    imputation: str | None
    ambient: str | None
    remove_doublets: bool
    doublet_methods: str | None
    condition_key: str | None
    deg_engine: str | None
    marker_method: str | None
    deg_cross_validate: str | bool | None
    qc_method: str | None
    report_lang: str
    phase: str  # qc | downstream | interpret
    metadata: dict[str, Any]
    plan: dict[str, Any]
    qc_strategy: dict[str, Any]
    annotation_plan: dict[str, Any]
    cluster_deg_plan: dict[str, Any]
    interpretation_plan: dict[str, Any]
    code_qc: str
    code_downstream: str
    code_interpret: str
    code: str
    execution_qc: dict[str, Any]
    execution_downstream: dict[str, Any]
    execution_interpret: dict[str, Any]
    execution: dict[str, Any]
    review_qc: dict[str, Any]
    review_downstream: dict[str, Any]
    review_publication: dict[str, Any]
    review: dict[str, Any]
    retry_count_qc: int
    retry_count_downstream: int
    artifacts: dict[str, Any]
    logs: list[str]
    report: str
    skills_used: list[str]
    r_degraded: bool
    status: str
    thread_id: str
    analysis_memory: dict[str, Any]
    notebook: str
    viewer: str
    selection: dict[str, Any]
    snapshots: list[dict[str, Any]]
