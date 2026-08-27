from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from agents.artifacts import skills_fingerprint
from scagent.config import analysis_params
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
    payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "skills_fingerprint": skills_fingerprint(),
        "seed": analysis_params()["seed"],
        **extra,
    }
    path = workspace / "run_manifest.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
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


def write_and_maybe_run(
    code: str,
    *,
    workspace: Path,
    execute: bool,
    timeout: int = 600,
    filename: str = "analysis.py",
    extra_manifest: dict | None = None,
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
    }
    if not execute or not code:
        result["stderr"] = "未执行代码。已写入脚本与 run_manifest.json。"
        return result

    needs_scanpy = filename in {"qc_preprocess.py", "cluster_annotate.py"} or "import scanpy" in code
    if needs_scanpy:
        missing = check_packages()
        result["missing_packages"] = missing
        if missing:
            result["ok"] = False
            result["stderr"] = "缺少依赖: " + ", ".join(missing) + "。pip install -r requirements-analysis.txt"
            return result

    env = os.environ.copy()
    env.setdefault("PYTHONHASHSEED", str(analysis_params()["seed"]))
    log.info("execute %s timeout=%s", script, timeout)
    with timed(f"execute.{filename}", log):
        proc = subprocess.run(
            [sys.executable, str(script)],
            cwd=str(workspace),
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
    result["executed"] = True
    result["stdout"] = proc.stdout
    result["stderr"] = proc.stderr
    result["returncode"] = proc.returncode
    result["ok"] = proc.returncode == 0
    result["figures"] = [str(p) for p in sorted((workspace / "figures").glob("*")) if p.is_file()]
    if not result["ok"]:
        log.error("script failed returncode=%s", proc.returncode)
    return result
