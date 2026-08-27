from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def write_and_maybe_run(
    code: str,
    *,
    workspace: Path,
    execute: bool,
    timeout: int = 600,
) -> dict:
    workspace.mkdir(parents=True, exist_ok=True)
    (workspace / "figures").mkdir(exist_ok=True)
    script = workspace / "analysis.py"
    script.write_text(code, encoding="utf-8")
    result = {
        "ok": True,
        "script": str(script),
        "stdout": "",
        "stderr": "",
        "executed": False,
        "figures": [str(p) for p in sorted((workspace / "figures").glob("*"))],
    }
    if not execute:
        result["stderr"] = "未执行代码（analysis.execute_code=false）。已写入 workspace/analysis.py。"
        return result
    proc = subprocess.run(
        [sys.executable, str(script)],
        cwd=str(workspace),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    result["executed"] = True
    result["stdout"] = proc.stdout
    result["stderr"] = proc.stderr
    result["ok"] = proc.returncode == 0
    result["figures"] = [str(p) for p in sorted((workspace / "figures").glob("*"))]
    return result
