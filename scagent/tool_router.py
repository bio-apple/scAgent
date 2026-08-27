"""SciAgent Tool Router: R-first defaults with Python fallback."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

from scagent.config import REPO_ROOT, load_config
from scagent.logutil import get_logger

log = get_logger("tool_router")

R_DIR = Path(__file__).resolve().parent / "r"

# module -> (R tool id, Python fallback id)
_BUILTIN_DEFAULTS: dict[str, tuple[str, str | None]] = {
    "qc": ("seurat", "scanpy"),
    "normalize": ("seurat", "scanpy"),
    "integration": ("harmony_r", "harmony_py"),
    "annotation": ("azimuth", "scanpy"),
    "trajectory": ("monocle3", "scanpy"),
    "cellchat": ("cellchat", None),
    "spatial": ("giotto", "squidpy"),
    "deg": ("edger", "ttest"),
}

_R_PACKAGES: dict[str, list[str]] = {
    "seurat": ["Seurat"],
    "harmony_r": ["harmony", "Seurat"],
    "azimuth": ["Azimuth", "Seurat"],
    "monocle3": ["monocle3"],
    "cellchat": ["CellChat", "Seurat"],
    "giotto": ["Giotto"],
    "edger": ["edgeR"],
}

_R_SCRIPTS: dict[str, str] = {
    "qc": "pipeline_qc.R",
    "normalize": "pipeline_qc.R",
    "integration": "harmony_integrate.R",
    "annotation": "pipeline_annotate.R",
    "trajectory": "monocle3.R",
    "cellchat": "cellchat.R",
    "spatial": "giotto_spatial.R",
}


def router_cfg(cfg: dict | None = None) -> dict:
    cfg = cfg or load_config()
    tr = dict(cfg.get("tool_router") or {})
    policy = str(tr.get("policy") or "r_first").lower().strip()
    defaults = dict(_BUILTIN_DEFAULTS)
    for mod, tool in (tr.get("defaults") or {}).items():
        fb = (tr.get("python_fallback") or {}).get(mod)
        if fb is None and mod in defaults:
            fb = defaults[mod][1]
        defaults[str(mod)] = (str(tool), str(fb) if fb else None)
    return {"policy": policy, "defaults": defaults, "raw": tr}


def analysis_language(cfg: dict | None = None) -> str:
    cfg = cfg or load_config()
    return str((cfg.get("analysis") or {}).get("language") or "r_first").lower().strip()


def rscript_path() -> str | None:
    return shutil.which("Rscript")


def r_packages_available(packages: list[str]) -> bool:
    r = rscript_path()
    if not r:
        return False
    pkgs = ", ".join(repr(p) for p in packages)
    cmd = f"all(vapply(c({pkgs}), requireNamespace, logical(1), quietly=TRUE))"
    try:
        proc = subprocess.run([r, "-e", cmd], capture_output=True, text=True, timeout=30, check=False)
        return proc.returncode == 0 and "TRUE" in (proc.stdout or "")
    except Exception as exc:
        log.debug("R package probe failed: %s", exc)
        return False


def r_tool_ready(tool_id: str) -> bool:
    pkgs = _R_PACKAGES.get(tool_id)
    if not pkgs:
        return bool(rscript_path())
    return r_packages_available(pkgs)


def resolve_module(
    module: str,
    *,
    cfg: dict | None = None,
    meta: dict | None = None,
    plan: dict | None = None,
) -> dict[str, Any]:
    """Pick primary engine for a module under current policy."""
    rc = router_cfg(cfg)
    policy = rc["policy"]
    lang = analysis_language(cfg)
    if lang == "r":
        policy = "r_only"
    elif lang == "python":
        policy = "python_only"

    defaults = rc["defaults"]
    primary, fallback = defaults.get(module, ("scanpy", None))
    meta = meta or {}
    plan = plan or {}

    if module == "integration" and not (meta.get("need_batch_correction") or plan.get("integrator")):
        return {
            "module": module,
            "engine": "none",
            "tool": "none",
            "fallback_tool": fallback,
            "policy": policy,
            "reason": "single sample / no integration",
        }
    if module == "spatial" and not _needs_spatial(meta, plan):
        return {
            "module": module,
            "engine": "none",
            "tool": "none",
            "fallback_tool": fallback,
            "policy": policy,
            "reason": "no spatial task",
        }
    if module == "cellchat" and "cellchat" not in (plan.get("route") or []) and "cellchat" not in str(plan.get("objective") or "").lower():
        if not any(k in str(meta.get("user_query") or plan.get("objective") or "").lower() for k in ("cellchat", "配体", "ligand", "通讯")):
            return {
                "module": module,
                "engine": "none",
                "tool": "none",
                "fallback_tool": fallback,
                "policy": policy,
                "reason": "cellchat not requested",
            }

    if policy == "python_only":
        return _pick(fallback or primary, "python", primary, fallback, policy, "python_only policy")

    if policy == "r_only" or policy == "r_first":
        if r_tool_ready(primary):
            return _pick(primary, "r", primary, fallback, policy, f"R tool {primary} available")
        if policy == "r_only":
            return _pick(fallback or "scanpy", "python", primary, fallback, policy, f"R tool {primary} missing; r_only degrades")
        if fallback:
            return _pick(fallback, "python", primary, fallback, policy, f"R tool {primary} missing; Python fallback")
        return _pick(primary, "r", primary, fallback, policy, f"R tool {primary} forced (no fallback)")

    # python_first
    if fallback:
        return _pick(fallback, "python", primary, fallback, policy, "python_first policy")
    return _pick(primary, "r" if r_tool_ready(primary) else "python", primary, fallback, policy, "python_first")


def _pick(tool: str, engine: str, primary: str, fallback: str | None, policy: str, reason: str) -> dict[str, Any]:
    return {
        "module": "",
        "engine": engine,
        "tool": tool,
        "primary_tool": primary,
        "fallback_tool": fallback,
        "policy": policy,
        "reason": reason,
        "rscript": rscript_path(),
    }


def _needs_spatial(meta: dict, plan: dict) -> bool:
    route = plan.get("route") or []
    if "spatial" in route:
        return True
    if meta.get("platform") in {"visium", "spatial", "xenium", "merfish"}:
        return True
    q = str(meta.get("user_query") or plan.get("objective") or "").lower()
    return any(k in q for k in ("spatial", "空间", "giotto", "squidpy", "visium"))


def build_tool_route(
    meta: dict | None = None,
    plan: dict | None = None,
    *,
    cfg: dict | None = None,
    modules: tuple[str, ...] | None = None,
) -> dict[str, Any]:
    """Full routing table for planner / report."""
    mods = modules or ("qc", "normalize", "integration", "annotation", "trajectory", "cellchat", "spatial", "deg")
    cfg = cfg or load_config()
    routes: dict[str, dict] = {}
    for mod in mods:
        r = resolve_module(mod, cfg=cfg, meta=meta, plan=plan)
        r["module"] = mod
        routes[mod] = r
    summary = {
        "policy": router_cfg(cfg)["policy"],
        "language": analysis_language(cfg),
        "r_available": bool(rscript_path()),
        "system_prompt": "Always use R ecosystem first. Only invoke Python when R lacks the required functionality.",
        "routes": routes,
    }
    return summary


def r_script_for_module(module: str) -> Path | None:
    name = _R_SCRIPTS.get(module)
    if not name:
        return None
    path = R_DIR / name
    return path if path.is_file() else None


def run_r_phase(
    module: str,
    *,
    workspace: Path,
    args: list[str],
    timeout: int = 600,
) -> dict[str, Any]:
    """Run an R pipeline script; return execution dict."""
    r = rscript_path()
    script = r_script_for_module(module)
    if not r or not script:
        return {"ok": False, "executed": False, "reason": "no_rscript_or_script", "engine": "r"}
    cmd = [r, str(script), *args]
    log.info("tool_router R: %s", " ".join(cmd))
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(workspace),
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return {
            "ok": False,
            "executed": True,
            "engine": "r",
            "returncode": -1,
            "stderr": str(exc),
            "stdout": "",
        }
    ok = proc.returncode == 0
    metrics_path = workspace / "r_metrics.json"
    metrics = {}
    if metrics_path.is_file():
        try:
            metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        except Exception:
            metrics = {}
    return {
        "ok": ok,
        "executed": True,
        "engine": "r",
        "returncode": proc.returncode,
        "stdout": proc.stdout or "",
        "stderr": proc.stderr or "",
        "metrics": metrics,
        "script": str(script),
    }


def format_route_table(route: dict[str, Any], *, lang: str = "zh") -> str:
    zh = lang != "en"
    lines = [
        "**Tool Router（R 优先）**" if zh else "**Tool Router (R-first)**",
        "",
        f"- policy: `{route.get('policy')}` | Rscript: `{route.get('r_available')}`",
        "",
        "| 功能 | 选用 | engine | 原因 |" if zh else "| module | tool | engine | reason |",
        "|---|---|---|---|",
    ]
    for mod, r in (route.get("routes") or {}).items():
        if r.get("engine") == "none":
            continue
        lines.append(
            f"| {mod} | {r.get('tool')} | {r.get('engine')} | {r.get('reason', '')} |"
        )
    lines.append("")
    lines.append(route.get("system_prompt") or "")
    return "\n".join(lines)
