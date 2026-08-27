"""Analysis provenance memory: steps and parameters, not chat."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml

from scagent.config import analysis_params, load_config, resolve_path


def _sample_id(state: dict) -> str:
    meta = state.get("metadata") or {}
    path = state.get("data_path") or meta.get("data_path") or ""
    stem = Path(str(path)).stem if path else ""
    return stem or str(meta.get("tissue") or "unknown")


def _step(code: str | None, execution: dict | None, snapshot: str | None) -> dict[str, Any]:
    execution = execution or {}
    if not code:
        status = "pending"
    elif not execution.get("executed"):
        status = "planned"
    elif execution.get("ok"):
        status = "ok"
    else:
        status = "failed"
    out: dict[str, Any] = {"status": status}
    if snapshot:
        out["snapshot"] = snapshot
    if execution.get("jail"):
        out["jail"] = execution.get("jail")
    return out


def _first_snapshot(execution: dict | None, h5ad: str | None) -> str | None:
    snaps = list((execution or {}).get("snapshots") or [])
    if snaps:
        return snaps[0]
    return h5ad or None


def _resume_from(steps: dict[str, dict]) -> str | None:
    for name in ("qc", "downstream"):
        st = (steps.get(name) or {}).get("status")
        if st in {"failed", "pending"}:
            return name
    return None


def build_memory(state: dict) -> dict[str, Any]:
    """Structured provenance for the run. No chat, no code blobs."""
    meta = state.get("metadata") or {}
    plan = state.get("plan") or {}
    qc = state.get("qc_strategy") or {}
    ann = state.get("annotation_plan") or {}
    arts = state.get("artifacts") or {}
    mets = arts.get("metrics") or {}
    h5ads = arts.get("h5ads") or {}
    hard = qc.get("hard") or {}
    params = analysis_params()
    exe_qc = state.get("execution_qc") or {}
    exe_dn = state.get("execution_downstream") or {}
    code_dn = state.get("code_downstream") or ""
    layer = meta.get("expression_layer")
    if layer in {"log1p", "scaled"}:
        normalize = f"skipped ({layer})"
    else:
        normalize = "normalize_total+log1p"

    annotation: list[str] = []
    if plan.get("celltypist_model") or "celltypist" in code_dn.lower():
        annotation.append("CellTypist")
    if ann.get("dual_validation") or ("positive" in code_dn.lower() and "negative" in code_dn.lower()):
        annotation.append("Marker")
    if "ref2_label" in code_dn or ann.get("ref2"):
        annotation.append("ref2")
    if not annotation:
        annotation = ["Marker"]

    if plan.get("needs_pseudobulk") or meta.get("needs_pseudobulk"):
        deg = "pseudobulk+FDR"
    elif "rank_genes_groups" in code_dn.lower() or "wilcox" in code_dn.lower():
        deg = "wilcox (exploratory)"
    else:
        deg = None

    steps = {
        "qc": _step(state.get("code_qc"), exe_qc, _first_snapshot(exe_qc, h5ads.get("qc"))),
        "downstream": _step(state.get("code_downstream"), exe_dn, _first_snapshot(exe_dn, h5ads.get("processed"))),
    }
    tid = state.get("thread_id")
    return {
        "sample": _sample_id(state),
        "tissue": meta.get("tissue"),
        "thread_id": tid,
        "qc": {
            "method": qc.get("method") or "mad",
            "nmads": qc.get("nmads"),
            "mt": mets.get("pct_mt_cutoff", hard.get("pct_mt")),
            "umi": mets.get("umi_min", hard.get("n_genes_min")),
            "n_before": mets.get("n_before"),
            "n_after": mets.get("n_after"),
            "pct_removed": mets.get("pct_removed"),
            "doublets": bool(qc.get("doublets", True)),
            "ambient": qc.get("ambient") or "none",
        },
        "normalize": normalize,
        "integration": plan.get("integrator"),
        "annotation": annotation,
        "deg": deg,
        "params": {
            "n_pcs": params["n_pcs"],
            "n_hvg": params["n_hvg"],
            "n_neighbors": params["n_neighbors"],
            "seed": params["seed"],
        },
        "steps": steps,
        "resume_from": _resume_from(steps),
        "resume": (
            f'python -m scagent run --from-checkpoint --thread-id {tid}'
            if tid
            else "python -m scagent run --from-checkpoint"
        ),
    }


def dump_memory_yaml(memory: dict) -> str:
    return yaml.safe_dump(memory, allow_unicode=True, sort_keys=False, default_flow_style=False)


def persist_memory(state: dict, extra_dir: Path | None = None, *, cfg: dict | None = None) -> dict[str, Any]:
    memory = build_memory(state)
    text = dump_memory_yaml(memory)
    cfg = cfg or load_config()
    cache = resolve_path(cfg, "cache")
    cache.mkdir(parents=True, exist_ok=True)
    (cache / "memory.yaml").write_text(text, encoding="utf-8")
    (cache / "memory.json").write_text(json.dumps(memory, ensure_ascii=False, indent=2), encoding="utf-8")
    if extra_dir is not None:
        extra_dir = Path(extra_dir)
        extra_dir.mkdir(parents=True, exist_ok=True)
        (extra_dir / "memory.yaml").write_text(text, encoding="utf-8")
    return memory


def load_memory(cfg: dict | None = None) -> dict[str, Any] | None:
    cfg = cfg or load_config()
    for path in (resolve_path(cfg, "cache") / "memory.yaml", resolve_path(cfg, "outputs") / "memory.yaml"):
        if path.is_file():
            data = yaml.safe_load(path.read_text(encoding="utf-8"))
            if isinstance(data, dict):
                return data
    return None
