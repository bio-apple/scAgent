from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

from agents.artifacts import parse_metrics, skills_fingerprint
from sandbox.jail import isolated_env, kernel_env, resolve_network, run_jailed, sandbox_settings
from sandbox.policy import policy_violations
from scagent.compat import scagent_version
from scagent.config import analysis_params, load_config
from scagent.logutil import get_logger, timed

log = get_logger("executor")

REQUIRED_PACKAGES = ("scanpy", "anndata", "numpy")


def check_packages(names: tuple[str, ...] = REQUIRED_PACKAGES) -> list[str]:
    missing = []
    for name in names:
        try:
            __import__(name)
        except ImportError:
            missing.append(name)
    return missing


def write_manifest(workspace: Path, extra: dict) -> Path:
    from scagent.reproducibility import enrich_run_manifest

    seed = analysis_params()["seed"]
    payload_extra = {
        "python": sys.version,
        "skills_fingerprint": skills_fingerprint(),
        "seed": seed,
        **extra,
        "scagent_version": scagent_version(),
    }
    path = workspace / "run_manifest.json"
    if path.is_file():
        try:
            prev = json.loads(path.read_text(encoding="utf-8"))
            payload_extra = {**prev, **payload_extra}
        except json.JSONDecodeError:
            pass
    enrich_run_manifest(path, extra=payload_extra)
    return path


def _refresh_reproducible(workspace: Path, filename: str, code: str) -> None:
    if filename == "qc_preprocess.py":
        (workspace / "reproducible_script.py").write_text(code, encoding="utf-8")
        return
    if filename == "cluster_annotate.py":
        p1 = workspace / "qc_preprocess.py"
        prefix = p1.read_text(encoding="utf-8") if p1.exists() else ""
        (workspace / "reproducible_script.py").write_text(
            prefix + "\n\n# --- PHASE 2 ---\n\n" + code, encoding="utf-8"
        )
        return
    (workspace / "reproducible_script.py").write_text(code, encoding="utf-8")


def _cleanup_tmp(workspace: Path) -> None:
    tmp = workspace / "sandbox_tmp"
    if tmp.is_dir():
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(exist_ok=True)


def analysis_executor(cfg: dict | None) -> str:
    raw = str(((cfg or {}).get("analysis") or {}).get("executor") or "jupyter").lower().strip()
    if raw in {"subprocess", "jail", "sandbox", "rlimit"}:
        return "subprocess"
    return "jupyter"


def _snapshot_h5ads(
    workspace: Path,
    phase: str,
    *,
    thread_id: str | None = None,
    params: dict | None = None,
    cfg: dict | None = None,
) -> list[dict]:
    from scagent.snapshot import record_phase

    return record_phase(workspace, phase, thread_id=thread_id, params=params, cfg=cfg)


def write_and_maybe_run(
    code: str,
    *,
    workspace: Path,
    execute: bool,
    timeout: int = 600,
    filename: str = "analysis.py",
    extra_manifest: dict | None = None,
    cfg: dict | None = None,
) -> dict:
    workspace.mkdir(parents=True, exist_ok=True)
    (workspace / "figures").mkdir(exist_ok=True)
    script = workspace / filename
    script.write_text(code or "", encoding="utf-8")
    if code:
        _refresh_reproducible(workspace, filename, code)

    extra_manifest = extra_manifest or {}
    data_path = extra_manifest.get("data_path")
    digest = None
    if data_path and Path(str(data_path)).is_file():
        h = hashlib.sha256()
        with open(data_path, "rb") as f:
            h.update(f.read(1024 * 1024))
        digest = h.hexdigest()

    write_manifest(
        workspace,
        {
            "script": str(script),
            "executed": bool(execute),
            "data_path": data_path,
            "data_sha256_head": digest,
            **extra_manifest,
        },
    )

    result = {
        "ok": True,
        "script": str(script),
        "stdout": "",
        "stderr": "",
        "executed": False,
        "returncode": 0,
        "figures": [str(p) for p in sorted((workspace / "figures").glob("*")) if p.is_file()],
        "missing_packages": [],
        "jail": None,
        "metrics": {},
        "warnings": [],
    }
    if not execute or not code:
        result["stderr"] = "未执行代码。已写入脚本与 run_manifest.json。"
        return result

    cfg = cfg or load_config()
    sb = sandbox_settings(cfg)
    phase = str((extra_manifest or {}).get("phase") or "")
    sb["network"] = resolve_network(sb, phase=phase)
    if sb.get("static_policy", True) and sb.get("enabled", True):
        blocked = policy_violations(code)
        if blocked:
            result["ok"] = False
            result["returncode"] = 126
            result["stderr"] = "sandbox policy blocked: " + ", ".join(blocked)
            result["jail"] = "policy"
            log.error("%s", result["stderr"])
            return result

    from agents.code_schema import validate_script

    schema = validate_script(code, phase=phase or "downstream")
    result["schema"] = {"ok": schema.get("ok"), "issues": schema.get("issues"), "steps": schema.get("steps")}
    if not schema.get("ok"):
        result["ok"] = False
        result["returncode"] = 125
        result["stderr"] = "schema validation failed: " + "; ".join(schema.get("issues") or ["invalid script"])
        result["jail"] = "schema"
        log.error("%s", result["stderr"])
        return result

    needs_scanpy = filename in {"qc_preprocess.py", "cluster_annotate.py"} or "import scanpy" in code
    if needs_scanpy:
        missing = check_packages()
        result["missing_packages"] = missing
        if missing:
            result["ok"] = False
            result["stderr"] = "缺少依赖: " + ", ".join(missing) + "。pip install -r requirements-analysis.txt"
            return result

    seed = analysis_params(cfg)["seed"]
    use_jupyter = analysis_executor(cfg) == "jupyter"
    if use_jupyter:
        env = kernel_env(workspace, seed=seed)
        env.setdefault("PYTHONHASHSEED", str(seed))
        log.info("execute %s timeout=%s executor=jupyter (no OS jail)", script, timeout)
        from scagent.export_nb import execute_via_jupyter

        with timed(f"execute.{filename}", log):
            ran = execute_via_jupyter(
                code,
                workspace=workspace,
                timeout=timeout,
                filename=filename,
                env=env,
                script=script,
            )
        result["executed"] = True
        result["stdout"] = ran.get("stdout") or ""
        result["stderr"] = ran.get("stderr") or ""
        result["returncode"] = int(ran.get("returncode") or 0)
        result["ok"] = bool(ran.get("ok"))
        result["jail"] = ran.get("jail") or "jupyter"
        result["notebook"] = ran.get("notebook")
        proc_returncode = result["returncode"]
        jail = result["jail"]
    else:
        env = isolated_env(workspace, seed=seed) if sb.get("enabled", True) else os.environ.copy()
        env.setdefault("PYTHONHASHSEED", str(seed))
        isolation = "off" if not sb.get("enabled", True) else sb.get("isolation") or "auto"
        sb = {**sb, "isolation": isolation}
        log.info("execute %s timeout=%s isolation=%s", script, timeout, isolation)
        with timed(f"execute.{filename}", log):
            proc = run_jailed(
                [sys.executable, str(script)],
                workspace=workspace,
                env=env,
                timeout=timeout,
                settings=sb,
            )
            jail = getattr(proc, "jail", isolation)
        result["executed"] = True
        result["stdout"] = proc.stdout
        result["stderr"] = proc.stderr
        result["returncode"] = proc.returncode
        result["ok"] = proc.returncode == 0
        result["jail"] = jail
        proc_returncode = proc.returncode
    result["figures"] = [str(p) for p in sorted((workspace / "figures").glob("*")) if p.is_file()]
    metrics, warns = parse_metrics(result["stdout"], result["stderr"])
    result["metrics"] = metrics
    result["warnings"] = warns
    if result.get("executed"):
        try:
            from scagent.reproducibility import enrich_run_manifest

            enrich_run_manifest(
                workspace / "run_manifest.json",
                extra={
                    "phase": phase or filename,
                    "metrics": metrics,
                    "data_path": extra_manifest.get("data_path"),
                    "executed": True,
                    "returncode": result.get("returncode"),
                },
            )
        except Exception as exc:
            log.warning("run_manifest provenance enrich failed: %s", exc)
    if result["ok"] and sb.get("enabled", True):
        from scagent.config import analysis_params as _aparams

        entries = _snapshot_h5ads(
            workspace,
            phase or filename,
            thread_id=str(extra_manifest.get("thread_id") or ""),
            params=_aparams(cfg),
            cfg=cfg,
        )
        result["snapshot_manifests"] = entries
        result["snapshots"] = [e.get("path") or e.get("obs") for e in entries if e]
        if sb.get("cleanup_tmp", True):
            _cleanup_tmp(workspace)
    if not result["ok"]:
        log.error("script failed returncode=%s jail=%s", proc_returncode, jail)
    return result
