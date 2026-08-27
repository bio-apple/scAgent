"""Run analysis phases via Tool Router (R-first, Python fallback)."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from typing import Any

from scagent.config import load_config
from scagent.logutil import get_logger
from scagent.tool_router import build_tool_route, resolve_module, run_r_phase

log = get_logger("phase_runner")


def _route(meta: dict, plan: dict, module: str) -> dict[str, Any]:
    return resolve_module(module, cfg=load_config(), meta=meta, plan=plan)


def maybe_run_r_qc(
    *,
    data_path: str,
    workspace: Path,
    meta: dict | None = None,
    plan: dict | None = None,
    qc: dict | None = None,
    timeout: int = 600,
) -> bool:
    """Try Seurat QC via Rscript. Return True if R path completed (caller should skip Python QC)."""
    if os.environ.get("SCAGENT_FORCE_PYTHON") == "1":
        return False
    meta = meta or {}
    plan = plan or {}
    qc = qc or {}
    route = _route(meta, plan, "qc")
    if route.get("engine") != "r":
        print("SCAGENT_TOOL_ROUTER qc engine=" + str(route.get("engine")) + " tool=" + str(route.get("tool")))
        return False
    ws = Path(workspace)
    ws.mkdir(parents=True, exist_ok=True)
    out = ws / "adata_qc.h5ad"
    sample_key = str(meta.get("sample_key") or "sample")
    nmads = int(qc.get("nmads") or 5)
    res = run_r_phase(
        "qc",
        workspace=ws,
        args=[str(data_path), str(out), sample_key, str(nmads)],
        timeout=timeout,
    )
    if not res.get("ok") or not out.is_file():
        tail = (res.get("stderr") or res.get("stdout") or "")[-500:]
        print("SCAGENT_WARN: R QC failed (" + tail + "); falling back to Scanpy")
        return False
    metrics = dict(res.get("metrics") or {})
    metrics["tool_router"] = route
    metrics["phase"] = "qc"
    print("SCAGENT_METRICS:" + json.dumps(metrics))
    (ws / "qc_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    cache = ws / ".cache"
    cache.mkdir(exist_ok=True)
    shutil.copy2(out, cache / "adata_qc.h5ad")
    print("SCAGENT_R_PHASE_COMPLETE qc")
    return True


def maybe_run_r_downstream(
    *,
    workspace: Path,
    meta: dict | None = None,
    plan: dict | None = None,
    timeout: int = 600,
) -> bool:
    """Try R annotate (+ optional harmony/cellchat). Return True if full R downstream succeeded."""
    if os.environ.get("SCAGENT_FORCE_PYTHON") == "1":
        return False
    meta = meta or {}
    plan = plan or {}
    ws = Path(workspace)
    inp = ws / "adata_qc.h5ad"
    if not inp.is_file():
        inp = ws / ".cache" / "adata_qc.h5ad"
    if not inp.is_file():
        return False

    working = ws / "adata_work.h5ad"
    shutil.copy2(inp, working)

    integ = _route(meta, plan, "integration")
    if integ.get("engine") == "r" and meta.get("need_batch_correction"):
        batch_key = str(meta.get("sample_key") or plan.get("sample_key") or "sample")
        out_i = ws / "adata_integrated.h5ad"
        res_i = run_r_phase(
            "integration",
            workspace=ws,
            args=[str(working), str(out_i), batch_key],
            timeout=timeout,
        )
        if res_i.get("ok") and out_i.is_file():
            working = out_i

    ann = _route(meta, plan, "annotation")
    if ann.get("engine") != "r":
        print("SCAGENT_TOOL_ROUTER annotate engine=" + str(ann.get("engine")))
        return False

    tissue = str(meta.get("tissue") or "default")
    resolution = plan.get("resolution") or 0.6
    out = ws / "adata_annotated.h5ad"
    res = run_r_phase(
        "annotation",
        workspace=ws,
        args=[str(working), str(out), tissue, str(resolution)],
        timeout=timeout,
    )
    if not res.get("ok") or not out.is_file():
        tail = (res.get("stderr") or res.get("stdout") or "")[-500:]
        print("SCAGENT_WARN: R annotation failed (" + tail + "); falling back to Scanpy")
        return False

    cc = _route(meta, plan, "cellchat")
    if cc.get("engine") == "r":
        run_r_phase("cellchat", workspace=ws, args=[str(out), "cellchat_out", "scagent_annotation"], timeout=timeout)

    sp = _route(meta, plan, "spatial")
    if sp.get("engine") == "r":
        run_r_phase("spatial", workspace=ws, args=[str(out), "spatial_out"], timeout=timeout)

    metrics = dict(res.get("metrics") or {})
    metrics["tool_router"] = {"annotation": ann, "integration": integ}
    metrics["phase"] = "downstream"
    print("SCAGENT_METRICS:" + json.dumps(metrics))
    (ws / "downstream_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    shutil.copy2(out, ws / "adata_qc.h5ad")
    print("SCAGENT_R_PHASE_COMPLETE downstream")
    return True


def route_for_state(state: dict) -> dict[str, Any]:
    return build_tool_route(state.get("metadata"), state.get("plan"))
